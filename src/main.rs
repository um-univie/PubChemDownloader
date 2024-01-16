use serde::Deserialize;
use std::io::prelude::*;
use tokio::time::{sleep, Duration};
#[derive(Deserialize, Debug)]
struct PubChemResponse {
    #[serde(rename = "IdentifierList")]
    record: Record,
}
#[derive(Deserialize, Debug)]
struct Record {
    #[serde(rename = "CID")]
    cid: Vec<u64>,
}

const MAX_CONCURRENT_REQUESTS: usize = 10;
const MAX_RETRIES: usize = 5;
const RETRY_DELAY: Duration = Duration::from_secs(30);

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = std::fs::read_to_string("molecules.txt")?;
    let smiles_vec: Vec<&str> = smiles.split('\n').collect();
    let molecules = std::fs::read_to_string("molecules.txt")?;
    let molecules_vec: Vec<&str> = molecules.split('\n').collect();
    let total_smiles = smiles_vec.len();
    let mut cids: Vec<Result<(String, u64), Box<dyn std::error::Error>>> = vec![];
    let mut smiles_list = vec![];

    let mut cids = vec![];
    let mut handles = vec![];
    for (index, smile) in smiles_vec.into_iter().enumerate() {
        if index % 5 == 0 {
            println!("Progress: {}%", 100 * index / total_smiles);
            println!(
                "Estimated time left: {} minutes",
                (total_smiles - index) * 5 / 60
            );
            cids.extend(futures::future::join_all(handles).await);
            handles = vec![];
        }
        handles.push(get_cid_by_name(smile));
    }

    let mut file = std::fs::File::create("cids.txt")?;

    for cid in cids.into_iter().flatten() {
        writeln!(file, "{},{}", cid.0, cid.1)?;
    }

    let mut cids = vec![];

    let file = std::fs::File::open("cids.txt")?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let mut line = line.split(',');
        let smiles = line.next().unwrap();
        let cid = line.next().unwrap().parse::<u64>().unwrap();
        cids.push(cid);
        smiles_list.push(smiles.to_string());
    }
    let mols = get_mols(cids, smiles_list);
    mols.await;
    println!("Molecules downloaded");
    Ok(())
}

async fn get_cid_by_smiles(smiles: &str) -> Result<(String, u64), Box<dyn std::error::Error>> {
    let url = format!(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/cids/JSON",
        smiles
    );
    let mut retry = 0;
    let retry_delay = Duration::from_secs(30);

    loop {
        match reqwest::get(&url).await {
            Ok(response) => {
                if response.status().is_success() {
                    let response_body = response.json::<PubChemResponse>().await?;
                    println!("{:?}", response_body);
                    println!("{}: {}", smiles, response_body.record.cid[0]);
                    return Ok((smiles.to_string(), response_body.record.cid[0]));
                } else if response.status().as_u16() == 429 {
                    // HTTP 429 Too Many Requests
                    println!("Rate limited for {}. Waiting to retry...", smiles);
                    sleep(retry_delay).await;
                } else {
                    return Err(format!("HTTP error {}: {}", response.status(), smiles).into());
                }
            }
            Err(e) => {
                println!("Error getting CID for {}: {}. Retrying...", smiles, e);
                sleep(retry_delay).await;
            }
        }

        retry += 1;
        if retry >= MAX_RETRIES {
            return Err(format!("Max retries reached for {}", smiles).into());
        }
    }
}

async fn get_cid_by_name(name: &str) -> Result<(String, u64), Box<dyn std::error::Error>> {
    let url = format!(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON",
        name
    );

    let mut retry = 0;
    let retry_delay = Duration::from_secs(30);

    println!("{:?}", url);
    loop {
        match reqwest::get(&url).await {
            Ok(response) => {
                if response.status().is_success() {
                    let response_body = response.json::<PubChemResponse>().await?;
                    println!("{:?}", response_body);
                    println!("{}: {}", name, response_body.record.cid[0]);
                    return Ok((name.to_string(), response_body.record.cid[0]));
                } else if response.status().as_u16() == 429 {
                    // HTTP 429 Too Many Requests
                    println!("Rate limited for {}. Waiting to retry...", name);
                    sleep(retry_delay).await;
                } else {
                    println!("{:?}", url);
                    println!("{:?}", response);
                    sleep(retry_delay).await;
                }
            }
            Err(e) => {
                println!("Error getting CID for {}: {}. Retrying...", name, e);
                sleep(retry_delay).await;
            }
        }

        retry += 1;
        if retry >= MAX_RETRIES {
            return Err(format!("Max retries reached for {}", name).into());
        }
    }
}

async fn get_mols(cids: Vec<u64>, molecule_names: Vec<String>) {
    for chunk in cids.chunks(MAX_CONCURRENT_REQUESTS) {
        if std::path::Path::new(&format!("sdf/{}.sdf", chunk[0])).exists() {
            continue;
        } else {
            std::fs::create_dir_all("sdf").unwrap();
        }
        let cid_string = chunk
            .iter()
            .map(|cid| cid.to_string())
            .collect::<Vec<String>>()
            .join(",");
        let url = format!(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{}/record/SDF?record_type=3d",
            cid_string
        );
        let handle = reqwest::get(&url).await.unwrap();
        println!("{}", url);
        let text = handle.text().await.unwrap();

        for sdf in text.split("$$$$") {
            let mut split_position = 0;
            if sdf.starts_with('\n') {
                split_position = 1;
            }
            let cid = match sdf.split('\n').nth(split_position) {
                Some(line) => {
                    println!("{line}");
                    line
                }
                None => continue,
            };
            if cid.is_empty() {
                continue;
            }
            let mut file =
                std::fs::File::create(format!("sdf/{cid}.sdf")).expect("Unable to create file");
            file.write_all(sdf.as_bytes()).unwrap();
        }
    }
}

async fn get_mol(cid: u64, molecule_name: &str) {
    let url = format!(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{}/record/SDF?record_type=3d",
        cid
    );
    let handle = reqwest::get(&url).await.unwrap();
    let text = handle.text().await.unwrap();
    let mut file =
        std::fs::File::create(format!("{cid}_{molecule_name}.sdf")).expect("Unable to create file");
    file.write_all(text.as_bytes()).unwrap();
}
