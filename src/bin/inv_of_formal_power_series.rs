use fps::*;
use std::io::*;
use std::collections::*;

fn main() {
  let mut input = String::new();
  stdin().read_to_string(&mut input).expect("An error occured while reading input");
  let mut words = input.split_ascii_whitespace().collect::<VecDeque<_>>();

  let N = words.pop_front().unwrap().parse::<usize>().unwrap();
  let a = (0 .. N).map(|_| words.pop_front().unwrap().parse::<i64>().unwrap() ).collect::<Vec<_>>();

  let f = FPS998244353::from(a);
  let g = f.inv_at(N);

  let ans = (0 .. N).map(|i| g[i].to_string() ).collect::<Vec<_>>().join(" ");
  println!("{}", ans);
}