use fps::FPS998244353 as FPS;
use acl_modint::ModInt998244353 as Mint;

fn main() {
  let mut f = FPS::from(vec![1, 1, 1, 1]) / FPS::from(vec![1, 1]);
  println!("{:?}", f);
}