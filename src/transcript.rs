use std::error::Error;
use std::collections::HashMap;
use noodles_gff::Record;
use crate::segment::Segment;

pub fn get_tid_from_rec(rec: gtf::Record) -> Option<String>{
    for i in 0..rec.attributes().len(){
        if rec.attributes()[i].key()=="transcript_id"{
            return Some(String::from(rec.attributes()[i].value()))
        }
    }
    None
}

#[derive(Debug)]
pub(crate) struct Transcript {
    records: Vec<Record>,
    tid: String,
}

impl Default for Transcript {
    fn default() -> Self {
        Transcript {
            records: vec![],
            tid: String,
        }
    }
}

impl Transcript {
    pub fn new(record: gtf::Record) -> Result<Self, Box<dyn Error>> {
        let tid = match get_tid_from_rec(record) {
            Some(tid) => tid,
            None => Err("transcript_id not fou"),
        };
        Ok(Self {
            records: vec![record],
            tid: tid,
        })
    }

    pub fn add(&mut self, record: gtf::Record) -> Result<(), Box<dyn Error>> {
        self.records.push(record);
        Ok(())
    }

    pub fn is_empty(&self) -> bool{
        self.records.is_empty()
    }
}