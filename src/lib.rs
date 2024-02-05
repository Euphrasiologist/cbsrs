/*!
The `cbsrs` crate implements the Circular Binary Segmentation algorithm (1)
based on code from (2).

It exposes an extremely simple API, as a trait which operates only on
`Vec<usize>`.

```no_run
use cbsrs::CBS;

let steps: Vec<usize> = vec![1, 1, 1, 3, 3, 2, 1, 2, 3, 300, 310, 321, 310, 299];

let shuffles = 1000;
let p = 0.05;

println!("{:?}", steps.cbs(shuffles, p));
```

This implementation omits the 'validation' algorithm seen in other implementations.

1. Olshen, Adam B., et al. "Circular binary segmentation for the analysis of array‐based DNA copy number data." Biostatistics 5.4 (2004): 557-572.
2. https://github.com/jeremy9959/cbs/blob/master/cbs.py
*/

use rand::seq::SliceRandom;
use rand::thread_rng;
use std::ops::Add;
use std::{error::Error as StdError, fmt, result};

pub type Result<T> = result::Result<T, Error>;

#[derive(Debug)]
pub struct Error(Box<ErrorKind>);

impl Error {
    /// A crate private constructor for `Error`.
    pub(crate) fn new(kind: ErrorKind) -> Error {
        Error(Box::new(kind))
    }

    /// Return the specific type of this error.
    pub fn kind(&self) -> &ErrorKind {
        &self.0
    }

    /// Unwrap this error into its underlying type.
    pub fn into_kind(self) -> ErrorKind {
        *self.0
    }
}

/// The kind of error possible, in this case
/// we only have one possibility.
#[derive(Debug)]
pub enum ErrorKind {
    /// Return the CBS error.
    CBS(String),
}

impl StdError for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &*self.0 {
            ErrorKind::CBS(err) => write!(f, "CBS failed - {}", err),
        }
    }
}

/// The `CBS` trait provides a single method, `cbs`, which
/// executes the circular binary segmentation algorithm.
pub trait CBS {
    /// Returns the circular binary segmentation output from the
    /// data. The output format is a vector of intervals (start, end)
    /// which slices up the data into distinct segments. User must specify a number of shuffles
    /// (1000 recommended), and a significance level (0.05 recommended).
    fn cbs(&self, shuffles: usize, significance_level: f64) -> Result<Vec<(usize, usize)>>;
}

impl CBS for Vec<usize> {
    /// The CBS implementation for a vector of type `usize`.
    fn cbs(&self, shuffles: usize, significance_level: f64) -> Result<Vec<(usize, usize)>> {
        let l = segment(self, shuffles, significance_level)?;

        Ok(l)
    }
}

fn cbs_stat(x: &[usize]) -> Result<(f64, usize, usize)> {
    if x.is_empty() {
        return Ok((0.0, 0, 0));
    }
    let length = x.len();
    let sum = x.iter().sum::<usize>();

    let mean = sum as f64 / length as f64;

    let x0: Vec<f64> = x.iter().map(|e| *e as f64 - mean).collect();

    let y = cumsum(&x0);

    let e0: usize = y
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .ok_or(Error::new(ErrorKind::CBS(
            "no maximum index in cbs_stat".into(),
        )))?;

    let e1: usize = y
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .ok_or(Error::new(ErrorKind::CBS(
            "no minimum index in cbs_stat".into(),
        )))?;

    let i0 = e0.min(e1);
    let i1 = e0.max(e1);

    let s0 = y[i0];
    let s1 = y[i1];

    // (s1-s0)**2*n/(i1-i0+1)/(n+1-i1+i0)
    Ok((
        (s1 - s0).powf(2.0) * length as f64
            / (i1 as f64 - i0 as f64 + 1.0)
            / (length as f64 - i1 as f64 + i0 as f64),
        i0,
        i1 + 1,
    ))
}

fn cbs_inner(x: &[usize], shuffles: usize, p: f64) -> Result<(bool, f64, usize, usize)> {
    // max_t, max_start, max_end = cbs_stat(x)
    let (max_t, mut max_start, mut max_end) = cbs_stat(x)?;
    // if max_end-max_start == len(x):
    if max_end - max_start == x.len() {
        //     return False, max_t, max_start, max_end
        return Ok((false, max_t, max_start, max_end));
    }

    // if max_start < 5:
    if max_start < 5 {
        //     max_start = 0
        max_start = 0;
    }
    // if len(x)-max_end < 5:
    if x.len() - max_end < 5 {
        //     max_end = len(x)
        max_end = x.len();
    }

    // thresh_count = 0
    let mut thresh_count = 0;
    // alpha = shuffles*p
    let alpha = shuffles as f64 * p;
    // xt = x.copy()
    let mut xt = x.to_vec();
    // for i in range(shuffles):
    for _ in 0..shuffles {
        //     np.random.shuffle(xt)
        xt.shuffle(&mut thread_rng());
        // threshold, s0, e0 = cbs_stat(xt)
        let (threshold, _, _) = cbs_stat(&xt)?;
        // if threshold >= max_t:
        if threshold >= max_t {
            // thresh_count += 1
            thresh_count += 1;
        }
        // if thresh_count > alpha:
        if thresh_count as f64 > alpha {
            // return False, max_t, max_start, max_end
            return Ok((false, max_t, max_start, max_end));
        }
    }
    // return True, max_t, max_start, max_end
    Ok((true, max_t, max_start, max_end))
}

fn rsegment(
    x: &[usize],
    start: usize,
    end: usize,
    l: &mut Vec<(usize, usize)>,
    shuffles: usize,
    p: f64,
) -> Result<Vec<(usize, usize)>> {
    let (threshold, _t, s, e) = cbs_inner(&x[start..end], shuffles, p)?;

    if !threshold || (e - s < 5) || (e - s == end - start) {
        l.push((start, end));
    } else {
        if s > 0 {
            rsegment(x, start, start + s, l, shuffles, p)?;
        }
        if e - s > 0 {
            rsegment(x, start + s, start + e, l, shuffles, p)?;
        }
        if start + e < end {
            rsegment(x, start + e, end, l, shuffles, p)?;
        }
    }

    let res = l.clone();

    Ok(res)
}

fn segment(x: &[usize], shuffles: usize, p: f64) -> Result<Vec<(usize, usize)>> {
    let start = 0;
    let end = x.len();
    let mut l = vec![];
    rsegment(x, start, end, &mut l, shuffles, p)?;

    Ok(l)
}

/// The cumulative sum of a generic vector, where elements
/// implement the add trait.
pub fn cumsum<T>(x: &[T]) -> Vec<T>
where
    T: Clone,
    for<'r> &'r T: Add<&'r T, Output = T>,
{
    let mut y = Vec::with_capacity(x.len());

    if !x.is_empty() {
        y.push(x[0].clone());

        for i in 1..x.len() {
            y.push(&y[i - 1] + &x[i]);
        }
    }

    y
}

#[cfg(test)]
mod tests {

    use super::CBS;

    #[test]
    fn test() {
        let steps: Vec<usize> = vec![1, 1, 1, 3, 3, 2, 1, 2, 3, 300, 310, 321, 310, 299];

        let shuffles = 1000;
        let p = 0.05;

        let res = steps.cbs(shuffles, p).unwrap();

        assert_eq!(res[0], (0, 8));
    }
}
