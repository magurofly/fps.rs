use acl_convolution::*;
use std::ops::*;

pub use acl_modint::*;
pub type FPS998244353 = FPS<ModInt998244353>;
pub type FPS1000000007 = FPS<ModInt1000000007>;

#[derive(Clone)]
/// Formal Power Series
pub struct FPS<T> {
  /// `a[i]` is coefficient of `x^i`
  a: Vec<T>,
}

impl<T: ModIntBase> FPS<T> {
  /// Returns `0`
  pub fn new() -> Self {
    Self { a: vec![] }
  }

  pub fn is_zero(&self) -> bool {
    self.len() == 0 || self.a.iter().all(|&x| x == T::from(0) )
  }

  pub fn constant(x: impl Into<T>) -> Self {
    Self { a: vec![x.into()] }
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

  pub fn differential(&self) -> Self {
    Self::from((1 .. self.len()).map(|i| self[i] * T::from(i) ))
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
    assert!(self.at(0) != StaticModInt::from(0));
    let mut ret = Self::new();
    ret[0] = self[0].inv();
    for i in 0 .. n.next_power_of_two().trailing_zeros() {
      ret = &ret + &ret - &ret * &ret * self.pre(1 << (i + 1));
      ret.a.truncate(1 << (i + 1));
    }
    ret.a.truncate(n);
    ret
  }

  pub fn integral(&self) -> Self {
    Self::from(Some(StaticModInt::from(0)).into_iter().chain((0 .. self.len()).map(|i| self[i] / StaticModInt::from(i + 1) )))
  }

  pub fn log(&self) -> Self {
    self.log_n(self.len())
  }

  pub fn log_n(&self, n: usize) -> Self {
    assert!(self.at(0) == StaticModInt::from(1));
    let mut ret = self.differential() * self.inv_n(n);
    ret.a.truncate(n - 1);
    ret.integral()
  }

  pub fn exp(&self) -> Self {
    self.exp_n(self.len())
  }

  pub fn exp_n(&self, n: usize) -> Self {
    // バグってそう HELP!
    assert!(self.at(0) == StaticModInt::from(0));
    let mut ret = Self::constant(1);
    for i in 0 .. n.next_power_of_two().trailing_zeros() {
      ret = &ret * (self.pre(2 << i) + StaticModInt::from(1) - ret.log_n(2 << i));
      ret.a.truncate(2 << i);
    }
    ret.a.truncate(n);
    ret
  }

  pub fn pow(&self, e: usize) -> Self {
    self.pow_n(self.len(), e)
  }

  pub fn pow_n(&self, n: usize, e: usize) -> Self {
    if e == 0 {
      return Self::constant(1);
    }

    // f^e = exp(e * log(f)) を利用する
    // ただし、 log(f) をするには f[0] == 1 でないといけないので変形する
    if let Some(first_nonzero) = (0 .. self.len()).find(|&i| self[i] != StaticModInt::from(0) ) {
      let mut f = self >> first_nonzero;
      f *= f[0].inv();
      f = f.log_n(n);
      f *= StaticModInt::from(e);
      f = f.exp_n(n - first_nonzero);
      f *= f[0].pow(e as u64);
      f <<= first_nonzero * e;
      f
    } else {
      Self::new()
    }
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


macro_rules! impl_op_from_assign {
  ($op:ident : $Op:ident from $op_assign:ident : $OpAssign:ident for <$($Params:ident),*> ($Self:ty, $Other:ty) where $($where:tt)*) => {
    // op_assign(&mut self, other)
    impl<$($Params),*> $OpAssign<$Other> for $Self where $($where)* {
      fn $op_assign(&mut self, other: $Other) {
        self.$op_assign(&other);
      }
    }

    // op(self, &other)
    impl<$($Params),*> $Op<&$Other> for $Self where $($where)* {
      type Output = $Self;
      fn $op(mut self, other: &$Other) -> Self::Output {
        self.$op_assign(other);
        self
      }
    }

    // op(self, other)
    impl<$($Params),*> $Op<$Other> for $Self where $($where)* {
      type Output = $Self;
      fn $op(self, other: $Other) -> Self::Output {
        self.$op(&other)
      }
    }

    // op(&self, &other)
    impl<$($Params),*> $Op<&$Other> for &$Self where $($where)* {
      type Output = $Self;
      fn $op(self, other: &$Other) -> Self::Output {
        self.clone().$op(other)
      }
    }

    // op(&self, other)
    impl<$($Params),*> $Op<$Other> for &$Self where $($where)* {
      type Output = $Self;
      fn $op(self, other: $Other) -> Self::Output {
        self.$op(&other)
      }
    }
  }
}

impl<T: ModIntBase> AddAssign<&Self> for FPS<T> {
  fn add_assign(&mut self, other: &Self) {
    self.reserve(other.len());
    for i in 0 .. other.len() {
      self.a[i] += other[i];
    }
  }
}
impl_op_from_assign!(add:Add from add_assign:AddAssign for <T> (FPS<T>, FPS<T>) where T: ModIntBase);

impl<T: ModIntBase> AddAssign<&T> for FPS<T> {
  fn add_assign(&mut self, other: &T) {
    self[0] += *other;
  }
}
impl_op_from_assign!(add:Add from add_assign:AddAssign for <T> (FPS<T>, T) where T: ModIntBase);

impl<T: ModIntBase> SubAssign<&Self> for FPS<T> {
  fn sub_assign(&mut self, other: &Self) {
    self.reserve(other.len());
    for i in 0 .. other.len() {
      self.a[i] -= other[i];
    }
  }
}
impl_op_from_assign!(sub:Sub from sub_assign:SubAssign for <T> (FPS<T>, FPS<T>) where T: ModIntBase);

impl<T: ModIntBase> SubAssign<&T> for FPS<T> {
  fn sub_assign(&mut self, other: &T) {
    self[0] -= *other;
  }
}
impl_op_from_assign!(sub:Sub from sub_assign:SubAssign for <T> (FPS<T>, T) where T: ModIntBase);

impl<M: Modulus> MulAssign<&Self> for FPS<StaticModInt<M>> {
  fn mul_assign(&mut self, other: &Self) {
    self.shrink();
    self.a = convolution(&self.a, &other.a);
  }
}
impl_op_from_assign!(mul:Mul from mul_assign:MulAssign for <M> (FPS<StaticModInt<M>>, FPS<StaticModInt<M>>) where M: Modulus);

impl<T: ModIntBase> MulAssign<&T> for FPS<T> {
  fn mul_assign(&mut self, other: &T) {
    for i in 0 .. self.len() {
      self[i] *= *other;
    }
  }
}
impl_op_from_assign!(mul:Mul from mul_assign:MulAssign for <T> (FPS<T>, T) where T: ModIntBase);

impl<M: Modulus> DivAssign<&Self> for FPS<StaticModInt<M>> {
  fn div_assign(&mut self, rhs: &Self) {
    *self *= rhs.inv_n(self.len().max(rhs.len()));
  }
}
impl_op_from_assign!(div:Div from div_assign:DivAssign for <M> (FPS<StaticModInt<M>>, FPS<StaticModInt<M>>) where M: Modulus);

impl<T: ModIntBase> ShrAssign<&usize> for FPS<T> {
  /// x^n 未満の項を消す
  fn shr_assign(&mut self, &n: &usize) {
    self.a.splice(0 .. n, None);
  }
}
impl_op_from_assign!(shr:Shr from shr_assign:ShrAssign for <T> (FPS<T>, usize) where T: ModIntBase);

impl<T: ModIntBase> ShlAssign<&usize> for FPS<T> {
  /// x^n をかける
  fn shl_assign(&mut self, &n: &usize) {
    self.a.splice(0 .. 0, (0 .. n).map(|_| T::from(0) ));
  }
}
impl_op_from_assign!(shl:Shl from shl_assign:ShlAssign for <T> (FPS<T>, usize) where T: ModIntBase);
