use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;

use base64::prelude::*;

const TIDK_DATABASE: &str = "tidk_database.csv";

// Define the structure of the tidk database
#[derive(Serialize, Deserialize, Debug)]
pub struct TelomereRepeatRow {
    #[serde(rename = "Phylum")]
    pub phylum: String,
    #[serde(rename = "Order")]
    pub order: String,
    #[serde(rename = "Family")]
    pub family: String,
    #[serde(rename = "Species")]
    pub species: String,
    #[serde(rename = "Telomeric repeat")]
    pub telomeric_repeat: String,
    #[serde(rename = "Notes")]
    pub notes: String,
    #[serde(rename = "Ref")]
    pub reference: String,
}

// Function to get the database path
pub fn get_database_path() -> Result<PathBuf> {
    let base_dir =
        dirs::data_dir().context("Could not determine the base directory for application data")?;

    // Create path for the application directory within the data dir
    let app_dir = base_dir.join("tidk");
    // Ensure that the directory exists
    fs::create_dir_all(&app_dir).context("Failed to create application data directory")?;

    Ok(app_dir.join(TIDK_DATABASE)) // Path to the dataset file
}

// Fetch the data from the github repository
// return the sha and the data
pub fn fetch_and_save_data() -> Result<()> {
    eprintln!("tidk build: fetching and saving data from the remote repository.");
    // Define the URL for the dataset
    // here we take from the actual repo
    let url =
        "https://api.github.com/repos/tolkit/telomeric-identifier/contents/clades/curated.csv";

    // Fetch the data
    let client = reqwest::blocking::Client::new();
    let mut headers = reqwest::header::HeaderMap::new();

    headers.insert("authorization", "<authorization>".parse().unwrap());
    headers.insert("user-agent", "CUSTOM_NAME/1.0".parse().unwrap());

    let response = client.get(url).headers(headers).send()?;
    let response_json = response.json::<serde_json::Value>()?;

    // FIXME: use this sha later on to do better version control?
    let _sha = response_json["sha"]
        .as_str()
        .context("No sha found in response")?;
    let mut content = response_json["content"]
        .as_str()
        .context("No content found in response")?
        .to_string();

    // whitespaces were messing up base64 decoding
    content.retain(|c| !c.is_whitespace());

    // decode the content. I am assuming it's standard base64 encoding.
    let data_u8 = BASE64_STANDARD.decode(content)?;
    let data = String::from_utf8(data_u8)?;

    // write the data to file - the name of which is the SHA commit
    let database_path = get_database_path()?;
    fs::write(&database_path, data)?;

    eprintln!("Data fetched and saved to: {}", database_path.display());

    Ok(())
}
