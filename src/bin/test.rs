use fps::FPS998244353 as FPS;
use fps::ModInt998244353 as Mint;

fn main() {
  let f = FPS::from(vec![0, 2, 1]);

  println!("{:?}", f);
  println!("f^3 = {:?}", f.pow_n(10, 3));
}