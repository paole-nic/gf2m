use core::ops::{Add, AddAssign, Mul, MulAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub struct Gf<const M: u32, const POLY: u64>(u64);

impl<const M: u32, const POLY: u64> Add for Gf<M, POLY> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Gf(self.0 ^ other.0)
    }
}

impl<const M: u32, const POLY: u64> AddAssign for Gf<M, POLY> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<const M: u32, const POLY: u64> Mul for Gf<M, POLY> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let masque = (1u64 << M) - 1;
        let mut a = self.0;
        let mut b = other.0;
        let mut result: u64 = 0;

        while b != 0 {
            if b & 1 != 0 {
                result ^= a;
            }
            let deb = a & (1u64 << (M - 1)) != 0;
            a = (a << 1) & masque;
            if deb {
                a ^= POLY;
            }
            b >>= 1;
        }
        Gf(result)
    }
}

impl<const M: u32, const POLY: u64> MulAssign for Gf<M, POLY> {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl<const M: u32, const POLY: u64> Gf<M, POLY> {
    pub fn new(v: u64) -> Self {
        Gf(v & ((1u64 << M) - 1))
    }

    pub fn square(self) -> Self {
        let mut etale: u128 = 0;
        for i in 0..M {
            if self.0 & (1u64 << i) != 0 {
                etale |= 1u128 << (2 * i);
            }
        }
        Gf(Self::reduire(etale))
    }

    fn reduire(mut v: u128) -> u64 {
        for i in (M..(2 * M - 1)).rev() {
            if v & (1u128 << i) != 0 {
                v ^= 1u128 << i;
                v ^= (POLY as u128) << (i - M);
            }
        }
        v as u64
    }

    pub fn inv(self) -> Self {
        assert!(self.0 != 0, "0");
        let mut a = self;
        let mut result = Self::new(1);

        for _ in 0..(M - 1) {
            a = a.square();
            result = result * a;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type G8 = Gf<8, 0b0001_1011>;

    // Add

    #[test]
    fn add_neutre() {
        for x in 0..256 {
            let a = G8::new(x);
            assert_eq!(a + G8::new(0), a);
        }
    }

    #[test]
    fn add_bin() {
        for x in 0..256 {
            let a = G8::new(x);
            assert_eq!(a + a, G8::new(0));
        }
    }

    #[test]
    fn add_com() {
        for x in 0..256 {
            for y in 0..256 {
                assert_eq!(G8::new(x) + G8::new(y), G8::new(y) + G8::new(x));
            }
        }
    }

    //mul

    #[test]
    fn mul_zero() {
        for x in 0..256 {
            let a = G8::new(x);
            assert_eq!(a * G8::new(0), G8::new(0));
        }
    }

    #[test]
    fn mul_neutre() {
        for x in 0..256 {
            let a = G8::new(x);
            assert_eq!(a * G8::new(1), a);
            assert_eq!(G8::new(1) * a, a);
        }
    }

    #[test]
    fn mul_com() {
        for x in 0..256 {
            for y in 0..256 {
                assert_eq!(G8::new(x) * G8::new(y), G8::new(y) * G8::new(x));
            }
        }
    }

    #[test]
    fn div_zero() {
        for x in 1..256 {
            for y in 1..256 {
                assert_ne!(G8::new(x) * G8::new(y), G8::new(0));
            }
        }
    }

    #[test]
    fn dist() {
        for x in 0..256 {
            for y in 0..256 {
                for z in 0..256 {
                    assert_eq!(
                        G8::new(z) * (G8::new(x) + G8::new(y)),
                        G8::new(z) * G8::new(x) + G8::new(z) * G8::new(y)
                    );
                }
            }
        }
    }

    #[test]
    fn carre() {
        for x in 0..256 {
            let a = G8::new(x);
            assert_eq!(a.square(), a * a);
        }
    }

    #[test]
    fn inverse() {
        for x in 1..256 {
            let a = G8::new(x);
            assert_eq!(a * a.inv(), G8::new(1));
        }
    }
}
