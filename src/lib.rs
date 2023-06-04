use acl_convolution::*;
use std::ops::*;

pub use acl_modint::*;
pub type FPS998244353 = FPS<ModInt998244353>;
pub type FPS1000000007 = FPS<ModInt1000000007>;

#[derive(Clone)]
/// Formal Power Series
pub struct FPS<T: ModIntBase> {
  /// `a[i]` is coefficient of `x^i`
  a: Vec<T>,
}

impl<T: ModIntBase> FPS<T> {
  /// Returns `0`
  pub fn new() -> Self {
    Self { a: vec![] }
  }

  /// Degree
  pub fn len(&self) -> usize {
    self.a.len()
  }

  /// Get first `n` coefficients
  pub fn pre(&self, n: usize) -> Self {
    Self::from((0 .. n).map(|i| self.at(i) ))
  }

  /// Get `n`-th coefficient
  pub fn at(&self, n: usize) -> T {
    if n < self.len() {
      self[n]
    } else {
      T::default()
    }
  }

  fn reserve(&mut self, i: usize) {
    if self.len() <= i {
      self.a.resize_with(i, T::default);
    }
  }

  fn shrink(&mut self) {
    let i = (0 .. self.len()).rev().find(|&i| self.a[i] != T::default() ).unwrap_or(0);
    self.a.truncate(i + 1);
  }
}

impl<M: Modulus> FPS<StaticModInt<M>> {
  /// Inverse
  pub fn inv(&self) -> Self {
    self.inv_n(self.len())
  }

  /// First `n` coefficients of inverse
  pub fn inv_n(&self, n: usize) -> Self {
    assert!(self.at(0) != Default::default());
    let mut ret = Self::new();
    ret[0] = self[0].inv();
    for i in 0 .. n.next_power_of_two().trailing_zeros() {
      ret = &ret + &ret - &ret * &ret * self.pre(1 << (i + 1));
      ret.a.truncate(1 << (i + 1));
    }
    ret.a.truncate(n);
    ret
  }
}

impl<T: ModIntBase> std::fmt::Debug for FPS<T> {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    f.write_fmt(format_args!("{:?}", self.a))
  }
}

impl<T: ModIntBase, U: Into<T>, I: IntoIterator<Item = U>> From<I> for FPS<T> {
  fn from(sequence: I) -> Self {
    Self { a: sequence.into_iter().map(U::into).collect() }
  }
}

impl<T: ModIntBase> Index<usize> for FPS<T> {
  type Output = T;
  fn index(&self, i: usize) -> &Self::Output {
    &self.a[i]
  }
}

impl<T: ModIntBase> IndexMut<usize> for FPS<T> {
  fn index_mut(&mut self, i: usize) -> &mut Self::Output {
    self.reserve(i + 1);
    &mut self.a[i]
  }
}

macro_rules! define_op_auxiliary {
    ($assign:ident, $Assign:ident, $op:ident, $Op:ident, <$($type:ident),*> [$($t:tt)*] where $($where:tt)*) => {
        impl<$($type)*> $Assign<Self> for FPS<$($t)*> where $($where)* {
          fn $assign(&mut self, other: Self) {
            self.$assign(&other);
          }
        }

        impl<$($type)*> $Op<&FPS<$($t)*>> for &FPS<$($t)*> where $($where)* {
          type Output = FPS<$($t)*>;
          fn $op(self, other: &FPS<$($t)*>) -> Self::Output {
            let mut this = self.clone();
            this.$assign(other);
            this
          }
        }

        impl<$($type)*> $Op<FPS<$($t)*>> for &FPS<$($t)*> where $($where)* {
          type Output = FPS<$($t)*>;
          fn $op(self, other: FPS<$($t)*>) -> Self::Output {
            let mut this = self.clone();
            this.$assign(other);
            this
          }
        }

        impl<$($type)*> $Op<&FPS<$($t)*>> for FPS<$($t)*> where $($where)* {
          type Output = FPS<$($t)*>;
          fn $op(mut self, other: &Self) -> Self::Output {
            self.$assign(other);
            self
          }
        }

        impl<$($type)*> $Op<FPS<$($t)*>> for FPS<$($t)*> where $($where)* {
          type Output = FPS<$($t)*>;
          fn $op(mut self, other: Self) -> Self::Output {
            self.$assign(other);
            self
          }
        }
    };
}

impl<T: ModIntBase> AddAssign<&Self> for FPS<T> {
  fn add_assign(&mut self, other: &Self) {
    self.reserve(other.len());
    for i in 0 .. other.len() {
      self.a[i] += other[i];
    }
  }
}
define_op_auxiliary!(add_assign, AddAssign, add, Add, <T> [T] where T: ModIntBase);

impl<T: ModIntBase> SubAssign<&Self> for FPS<T> {
  fn sub_assign(&mut self, other: &Self) {
    self.reserve(other.len());
    for i in 0 .. other.len() {
      self.a[i] -= other[i];
    }
  }
}
define_op_auxiliary!(sub_assign, SubAssign, sub, Sub, <T> [T] where T: ModIntBase);

impl<M: Modulus> MulAssign<&Self> for FPS<StaticModInt<M>> {
  fn mul_assign(&mut self, rhs: &Self) {
    self.shrink();
    self.a = convolution(&self.a, &rhs.a);
  }
}
define_op_auxiliary!(mul_assign, MulAssign, mul, Mul, <M> [StaticModInt<M>] where M: Modulus);

impl<M: Modulus> DivAssign<&Self> for FPS<StaticModInt<M>> {
  fn div_assign(&mut self, rhs: &Self) {
    *self *= rhs.inv_n(self.len().max(rhs.len()));
  }
}
define_op_auxiliary!(div_assign, DivAssign, div, Div, <M> [StaticModInt<M>] where M: Modulus);

#[cfg(test)]
mod tests {
}
