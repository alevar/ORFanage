use std::fs::File;
use std::error::Error;
use std::io::{BufRead, BufReader, Read};

// should be able to take either gtf or a vector of gtf files
#[derive(Default)]
pub(crate) struct TReader {
    fnames: Vec<String>,
    readers: Vec<BufReader<File>>,
    current_txs : Vec<Option<u32>>,
}

impl TReader {
    pub fn new<A>(args: A) -> TReader
        where A: Into<TReader>
    {
        args.into()
    }

    pub fn add(&mut self, fname: &String) -> Result<(),Box<dyn Error>>{
        self.fnames.push(fname.clone());
        self.current_txs.push(None);
        let file = File::open(fname)?;
        let mut reader = BufReader::new(file);

        // load the first value from this file
        let mut line = String::new();
        let res = reader.read_line(&mut line);
        if res.is_ok() && !line.is_empty() {
            let val = line.trim().parse::<u32>().expect("Failed to parse integer");
            self.current_txs[self.fnames.len() - 1] = Some(val);
        }

        self.readers.push(reader);

        Ok(())
    }
}

impl Iterator for TReader {
    // needs to
    type Item = Vec<u32>;
    fn next(&mut self) -> Option<Self::Item>{
        // get minimum value
        let min_value = self.current_txs.iter().min();
        let min_value = match min_value {
            Some(Some(min)) => *min,
            Some(None) => {return None;},
            None      => { return None; }
        };

        let mut group:Vec<u32> = vec![];

        for i in 0..self.current_txs.len(){
            if self.current_txs[i]==Some(min_value){
                group.push(min_value);
                loop{ // get new data until the value is no longer == to min_value
                    let mut new_value = None;
                    let mut line = String::new();
                    let res = self.readers[i].read_line(&mut line);
                    if res.is_ok() {
                        line = line.trim().to_string();
                        if !line.is_empty() {
                            let value = line.parse::<u32>().expect("Failed to parse integer");
                            new_value = Some(value);
                        }
                    }
                    if new_value!=Some(min_value){
                        self.current_txs[i] = new_value;
                        break;
                    }
                    else{
                        group.push(min_value);
                    }
                }
            }
        }

        Some(group)
    }
}

impl From<&String> for TReader {
    fn from(fname: &String) -> Self {
        let mut treader = TReader::default();
        treader.add(fname).expect("Unable to add file}");
        treader
    }
}

impl From<&Vec<&String>> for TReader {
    fn from(fnames: &Vec<&String>) -> Self {
        let mut treader = TReader::default();
        for f in fnames{
            treader.add(f).expect("Unable to add file}");
        }
        treader
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}