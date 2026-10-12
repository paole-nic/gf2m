use crate::gf::Gf;

#[derive(Debug, Clone, PartialEq, Eq)]

pub struct Vecteur<const M: u32, const POLY: u64> {
    coords: Vec<Gf<M, POLY>>,
}

impl<const M: u32, const POLY: u64> Vecteur<M, POLY> {
    pub fn from_vect(coords: Vec<Gf<M, POLY>>) -> Self {
        Vecteur { coords }
    }

    pub fn len(&self) -> usize {
        self.coords.len()
    }

    pub fn get(&self, i: usize) -> Gf<M, POLY> {
        self.coords[i]
    }

    pub fn sum(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len(), "longueur invalide");

        let mut coords = Vec::new();
        for i in 0..self.len() {
            coords.push(self.coords[i] + other.coords[i]);
        }
        Self { coords }
    }

    pub fn scalaire(&self, k: Gf<M, POLY>) -> Self {
        let mut coords = Vec::new();
        for i in 0..self.len() {
            coords.push(k * self.coords[i]);
        }
        Self { coords }
    }

    pub fn dot(&self, other: &Self) -> Gf<M, POLY> {
        assert_eq!(self.len(), other.len(), "longueur invalide");
        let mut a = Gf::<M, POLY>::new(0);
        for i in 0..self.len() {
            a += self.coords[i] * other.coords[i];
        }
        a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type G8 = Gf<8, 0b0001_1011>;
    type V8 = Vecteur<8, 0b0001_1011>;

    #[test]
    fn long() {
        let a = V8::from_vect(vec![G8::new(0), G8::new(0), G8::new(0)]);
        assert_eq!(a.len(), 3);
    }

    #[test]
    fn vec_get() {
        let a = V8::from_vect(vec![G8::new(1), G8::new(2), G8::new(3)]);
        for i in 1..4 {
            assert_eq!(a.get(i - 1), G8::new(i as u64));
        }
    }

    #[test]
    fn som_vec() {
        let a = V8::from_vect(vec![G8::new(1), G8::new(2), G8::new(3)]);
        let somme = a.sum(&a);
        assert_eq!(
            somme,
            V8::from_vect(vec![G8::new(0), G8::new(0), G8::new(0)])
        );
    }

    #[test]
    fn scal_vec() {
        let a = V8::from_vect(vec![G8::new(1), G8::new(2), G8::new(3)]);
        let k = G8::new(1);
        assert_eq!(a.scalaire(k), a);
    }

    #[test]
    fn prod_scal() {
        let a = V8::from_vect(vec![G8::new(1), G8::new(2), G8::new(3)]);
        let scal = a.dot(&a);
        let calc = G8::new(1).square() + G8::new(2).square() + G8::new(3).square();
        assert_eq!(scal, calc);
    }
}
