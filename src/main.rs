
extern crate nalgebra;
extern crate log;

use nalgebra::{
  DMatrix,
  Matrix,
  Complex,
  Dynamic,
  VecStorage,
  U1
};

use std::time::{
  Instant
};

fn solve_companion_matrix(coeffs: Vec<i64>, offset: i64) -> Matrix<Complex<f64>, Dynamic, U1, VecStorage<Complex<f64>, Dynamic, U1>> {
  let mut matrix = DMatrix::zeros(coeffs.len(), coeffs.len());

  for (ith, coeff) in coeffs.iter().enumerate() {
    if ith > 0 {
      matrix[(ith, ith - 1)] = 1 as f64;
    }

    matrix[(ith, coeffs.len() - 1)] = (coeff - offset) as f64;
  }

  let decomp = matrix.schur();
  return decomp.complex_eigenvalues();
}

fn to_mixed_radix(bases: Vec<i64>, mut num: i64) -> Vec<i64> {
  let mut ith = bases.len() - 1;
  let mut acc:Vec<i64> = vec![0; bases.len()];

  loop {
    acc[ith] = num % bases[ith];
    num = num / bases[ith];

    if ith == 0 {
      break;
    }

    ith = ith - 1;
  }

  return acc;
}

fn main() {
  let mut solved = 0;
  let now = Instant::now();

  loop {
    let eigens = solve_companion_matrix(to_mixed_radix(vec![50, 50, 50, 50], 7), solved);

    if solved > 1000 && solved % 10000 == 0 {
    let seconds_elapsed = now.elapsed().as_secs();
      let per_second = (solved as u64) / seconds_elapsed;
      println!("{} per second", per_second);
    }

    solved = solved + 1;
  }
}
