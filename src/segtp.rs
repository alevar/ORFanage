use std::error::Error;

#[derive(PartialEq)]
#[derive(Debug)]
pub(crate) struct SEGTP {
    start: u32,
    end: u32,
}

impl Default for SEGTP {
    fn default() -> Self {
        SEGTP {
            start: u32::MAX,
            end: u32::MAX,
        }
    }
}

impl SEGTP {
    pub fn new(start:u32, end:u32) -> Result<Self, Box<dyn Error>>{
        if start > end{
            return Err("start > end".into());
        }
        Ok(Self {
            start,
            end,
        })
    }

    pub fn set_start(&mut self,start:u32) -> Result<(),Box<dyn Error>>{
        if start > self.end{
            return Err("start > end".into());
        }
        self.start = start;
        Ok(())
    }

    pub fn set_end(&mut self,end:u32) -> Result<(),Box<dyn Error>>{
        if self.start > end{
            return Err("start > end".into());
        }
        self.end = end;
        Ok(())
    }

    pub fn contains(&self,pos: u32) -> bool{
        return !(pos < self.start || pos > self.end);
    }

    pub fn slen(&self) -> u32{
        return (self.end + 1) - self.start;
    }

    pub fn intersect(&self, s: &SEGTP) -> Option<SEGTP>{
        let start = self.start.max(s.start);
        let end = self.end.min(s.end);

        if start <= end {
            Some(SEGTP::new(start, end).unwrap())
        } else {
            None
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_segtp() {
        let segtp = SEGTP::new(1, 5).unwrap();
        assert_eq!(segtp.start, 1);
        assert_eq!(segtp.end, 5);
    }

    #[test]
    fn test_create_invalid_segtp() {
        assert!(SEGTP::new(5, 1).is_err());
    }

    #[test]
    fn test_set_start() {
        let mut segtp = SEGTP::new(1, 5).unwrap();
        assert!(segtp.set_start(0).is_ok());
        assert!(segtp.set_start(6).is_err());
        assert_eq!(segtp.start, 0);
        assert!(segtp.set_start(2).is_ok());
        assert_eq!(segtp.start, 2);
    }

    #[test]
    fn test_set_end() {
        let mut segtp = SEGTP::new(1, 5).unwrap();
        assert!(segtp.set_end(0).is_err());
        assert!(segtp.set_end(6).is_ok());
        assert_eq!(segtp.end, 6);
        assert!(segtp.set_end(4).is_ok());
        assert_eq!(segtp.end, 4);
    }

    #[test]
    fn test_contains() {
        let segtp = SEGTP::new(1, 5).unwrap();
        assert!(segtp.contains(1));
        assert!(segtp.contains(5));
        assert!(segtp.contains(3));
        assert!(!segtp.contains(0));
        assert!(!segtp.contains(6));
    }

    #[test]
    fn test_slen() {
        let segtp = SEGTP::new(1, 5).unwrap();
        assert_eq!(segtp.slen(), 5);
    }

    #[test]
    fn test_intersect() {
        let segtp1 = SEGTP::new(1, 5).unwrap();
        let segtp2 = SEGTP::new(3, 7).unwrap();
        let segtp3 = SEGTP::new(6, 10).unwrap();
        let segtp4 = SEGTP::new(0, 3).unwrap();
        let segtp5 = SEGTP::new(7, 9).unwrap();

        assert_eq!(SEGTP::intersect(&segtp1, &segtp2), Some(SEGTP::new(3, 5).unwrap()));
        assert_eq!(segtp1.intersect(&segtp3), None);
        assert_eq!(segtp1.intersect(&segtp4), Some(SEGTP::new(1, 3).unwrap()));
        assert_eq!(segtp1.intersect(&segtp5), None);
    }

}