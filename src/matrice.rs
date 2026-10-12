use crate::gf::Gf;
use crate::vect::Vecteur;

#[derive(Debug, Clone, PartialEq, Eq)]

pub struct Matrice<const M: u32, const POLY: u64> {
    lignes: usize,
    colones: usize,
    donnees: Vec<Gf<M, POLY>>,
}

impl<const M: u32, const POLY: u64> Matrice<M, POLY> {
    pub fn zero(lignes: usize, colones: usize) -> Self {
        Matrice {
            lignes,
            colones,
            donnees: vec![Gf::new(0); lignes * colones],
        }
    }

    pub fn get(&self, lignes: usize, colones: usize) -> Gf<M, POLY> {
        self.donnees[lignes * self.colones + colones]
    }

    pub fn set(&mut self, lignes: usize, colones: usize, a: Gf<M, POLY>) {
        let index = lignes * self.colones + colones;
        self.donnees[index] = a;
    }

    pub fn identite(n: usize) -> Self {
        let mut donnees = vec![Gf::new(0); n * n];
        for i in 0..n {
            donnees[i * n + i] = Gf::new(1);
        }
        Matrice {
            lignes: n,
            colones: n,
            donnees,
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.lignes, other.lignes, "nbr de lignes incompatible");
        assert_eq!(self.colones, other.colones, "nbr de colones incompatible");
        let mut donnees = vec![Gf::new(0); self.lignes * self.colones];
        for i in 0..(self.lignes * self.colones) {
            donnees[i] = self.donnees[i] + other.donnees[i];
        }
        Matrice {
            lignes: self.lignes,
            colones: self.colones,
            donnees,
        }
    }

    pub fn t(&self) -> Self {
        let mut donnees = vec![Gf::new(0); self.lignes * self.colones];
        for i in 0..self.colones {
            for j in 0..self.lignes {
                donnees[i * self.lignes + j] = self.donnees[j * self.colones + i];
            }
        }

        Matrice {
            lignes: self.colones,
            colones: self.lignes,
            donnees,
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.colones, other.lignes, "dimensions incompatibles");
        let mut donnees = vec![Gf::new(0); self.lignes * other.colones];

        for i in 0..self.lignes {
            for j in 0..other.colones {
                let mut a = Gf::new(0);
                for h in 0..self.colones {
                    a += self.donnees[i * self.colones + h] * other.donnees[h * other.colones + j];
                }
                donnees[i * other.colones + j] = a;
            }
        }

        Matrice {
            lignes: self.lignes,
            colones: other.colones,
            donnees,
        }
    }

    pub fn mul_mv(&self, other: &Vecteur<M, POLY>) -> Vecteur<M, POLY> {
        let mut coords = Vec::with_capacity(self.lignes);
        for i in 0..self.lignes {
            let mut v = Gf::new(0);
            for j in 0..self.colones {
                v += self.get(i, j) * other.get(j);
            }
            coords.push(v);
        }
        Vecteur::from_vect(coords)
    }

    pub fn pivot(&self) -> (Self, usize) {
        let mut m = self.clone();
        let mut ligne = 0;

        for col in 0..m.colones {
            if ligne >= m.lignes {
                break;
            }

            //capture ligne pivot
            let mut p = None;
            for i in ligne..m.lignes {
                if m.get(i, col) != Gf::new(0) {
                    p = Some(i);
                    break;
                }
            }
            let p = match p {
                Some(i) => i,
                None => continue,
            };

            //copie
            if p != ligne {
                for c in 0..m.colones {
                    let tmp = m.get(ligne, c);
                    m.set(ligne, c, m.get(p, c));
                    m.set(p, c, tmp);
                }
            }

            // Normalisation
            let inv = m.get(ligne, col).inv();
            for c in col..m.colones {
                let v = m.get(ligne, c) * inv;
                m.set(ligne, c, v);
            }

            // somme
            for i in 0..m.lignes {
                if i == ligne {
                    continue;
                }
                let facteur = m.get(i, col);
                if facteur == Gf::new(0) {
                    continue;
                }
                for c in col..m.colones {
                    let v = m.get(i, c) + facteur * m.get(ligne, c);
                    m.set(i, c, v);
                }
            }
            ligne += 1;
        }
        (m, ligne)
    }

    pub fn inv(&self) -> Self {
        assert_eq!(self.lignes, self.colones, "pas carre");
        let mut calc = Self::zero(self.lignes, 2 * self.colones);
        let mut m = Self::zero(self.lignes, self.colones);
        for l in 0..self.lignes {
            for c in 0..self.colones {
                calc.donnees[l * 2 * self.colones + c] = self.donnees[l * self.colones + c];
            }
            calc.donnees[l * 2 * self.lignes + l + self.lignes] = Gf::new(1);
        }
        let (calc, r) = calc.pivot();
        assert_eq!(r, self.lignes, "pas inversible");
        for l in 0..self.lignes {
            for c in 0..self.colones {
                m.donnees[l * self.colones + c] =
                    calc.donnees[l * 2 * self.colones + c + self.colones];
            }
        }
        m
    }

    pub fn count_pivot(&self) -> Vec<usize> {
        let mut p = Vec::with_capacity(self.lignes.min(self.colones));
        let mut c = 0;

        for l in 0..self.lignes {
            while c < self.colones && self.get(l, c) == Gf::new(0) {
                c += 1;
            }
            if c == self.colones {
                break;
            }
            p.push(c);
            c += 1;
        }
        p
    }

    pub fn ker(&self) -> Self {
        let (reduit, rang) = self.pivot();
        let pivots = reduit.count_pivot();

        let libres: Vec<usize> = (0..self.colones).filter(|c| !pivots.contains(c)).collect();

        let mut ker = Self::zero(self.colones, self.colones - rang);

        for (k, &j) in libres.iter().enumerate() {
            ker.set(j, k, Gf::new(1));
            for i in 0..rang {
                ker.set(pivots[i], k, reduit.get(i, j));
            }
        }

        ker
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type G8 = Gf<8, 0b0001_1011>;
    type V8 = Vecteur<8, 0b0001_1011>;
    type M8 = Matrice<8, 0b0001_1011>;

    #[test]
    fn set_get_id() {
        let mut m = M8::identite(2);
        m.set(0, 1, G8::new(2));
        assert_eq!(m.get(0, 0), G8::new(1));
        assert_eq!(m.get(0, 1), G8::new(2));
        assert_eq!(m.get(1, 0), G8::new(0));
        assert_eq!(m.get(1, 1), G8::new(1));
    }

    #[test]
    fn zero_mat() {
        let m = M8::zero(3, 4);
        for i in 0..3 {
            for j in 0..4 {
                assert_eq!(m.get(i, j), G8::new(0));
            }
        }
    }

    #[test]
    fn transpose() {
        let mut m = M8::zero(2, 3);
        for i in 0..2 {
            for j in 0..3 {
                m.set(i, j, G8::new((i * 3 + j) as u64));
            }
        }
        let t = m.t();
        assert_eq!(t.lignes, 3);
        assert_eq!(t.colones, 2);
        for i in 0..2 {
            for j in 0..3 {
                assert_eq!(t.get(j, i), m.get(i, j));
            }
        }
    }

    #[test]
    fn id_neutre() {
        let mut m = M8::zero(3, 3);
        for i in 0..3 {
            for j in 0..3 {
                m.set(i, j, G8::new((i * 3 + j + 1) as u64));
            }
        }
        let id = M8::identite(3);
        assert_eq!(m.mul(&id), m);
        assert_eq!(id.mul(&m), m);
    }

    #[test]
    fn transpose_produit() {
        let mut a = M8::zero(2, 3);
        let mut b = M8::zero(3, 2);
        for i in 0..2 {
            for j in 0..3 {
                a.set(i, j, G8::new((i * 3 + j + 1) as u64));
            }
        }
        for i in 0..3 {
            for j in 0..2 {
                b.set(i, j, G8::new((i * 2 + j + 5) as u64));
            }
        }
        assert_eq!(a.mul(&b).t(), b.t().mul(&a.t()));
    }

    #[test]
    fn rang() {
        let a = M8::zero(4, 4);
        let (_, rang) = a.pivot();
        assert_eq!(rang, 0);
    }

    #[test]
    fn pivot_id() {
        let id = M8::identite(3);
        let (r, rang) = id.pivot();
        assert_eq!(r, id);
        assert_eq!(rang, 3);
    }

    #[test]
    fn rang_matrice_sing() {
        let mut m = M8::zero(3, 3);
        for j in 0..3 {
            m.set(0, j, G8::new((j + 1) as u64));
            m.set(1, j, G8::new((j + 5) as u64));
            m.set(2, j, G8::new((j + 1) as u64));
        }
        let (_, rang) = m.pivot();
        assert_eq!(rang, 2);
    }

    #[test]
    fn ker_mat() {
        let mut m = M8::zero(3, 6);
        for i in 0..3 {
            for j in 0..6 {
                m.set(i, j, G8::new((j * 3 + 6 * i) as u64));
            }
        }
        assert_eq!(m.mul(&m.ker()), M8::zero(3, 3));
    }

    #[test]
    fn inv_me() {
        let mut m = M8::zero(3, 3);
        for i in 0..3 {
            for j in 0..3 {
                m.set(i, j, G8::new((i * 5 + j * 3 + 1) as u64));
            }
        }
        assert_eq!(m.mul(&m.inv()), M8::identite(3));
        assert_eq!(m.inv().mul(&m), M8::identite(3));
    }

    #[test]
    fn matrice_vecteur_identite() {
        let v = V8::from_vect(vec![G8::new(3), G8::new(7), G8::new(11)]);
        assert_eq!(M8::identite(3).mul_mv(&v), v);
    }
}
