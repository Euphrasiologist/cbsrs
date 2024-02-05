# `cbsrs`

The `cbsrs` crate implements the Circular Binary Segmentation algorithm (1) based on code from (2).

It exposes an extremely simple API, as a trait which operates only on `Vec<usize>`.

```rust
use cbsrs::CBS;

let steps: Vec<usize> = vec![1, 1, 1, 3, 3, 2, 1, 2, 3, 300, 310, 321, 310, 299];

let shuffles = 1000;
let p = 0.05;

let res = steps.cbs(shuffles, p).unwrap();

for (start, end) in res.iter() {
  println!("{start}-{end}");
}
```

This implementation omits the 'validation' algorithm seen in other implementations.

1. Olshen, Adam B., et al. "Circular binary segmentation for the analysis of array‐based DNA copy number data." Biostatistics 5.4 (2004): 557-572.
2. https://github.com/jeremy9959/cbs/blob/master/cbs.py
