use rand::Rng;

const N: usize = 150;
const K: usize = 2;

pub fn init_1d(cplngs: &mut Vec<Vec<f32>>) {
    for i in 1..N {
        let ii = i * K;
        for j in 0..K {
            cplngs[ii + j][ii - K + j] = -1.0;
            cplngs[ii - K + j][ii + j] = -1.0;
        }
    }
    for j in 0..K {
        cplngs[j][K + j] = -1.0;
        cplngs[K + j][j] = -1.0;
    }
}

pub fn delta_energy(system: &Vec<i8>, cplngs: &Vec<Vec<f32>>, pos: usize, old: usize, new: usize) -> f32 {
    let mut en: f32 = 0.0;
    let mut pos_j: usize = 0;
    let pos_i: usize = pos * K;
    for aa_j in system.iter() {
        en += cplngs[pos_j + *aa_j as usize][pos_i + new] - cplngs[pos_j + *aa_j as usize][pos_i + old];
        // println!("{} {} {} {} {}", pos_j, aa_j, cplngs[pos_j + *aa_j as usize][pos_i + new], cplngs[pos_j + *aa_j as usize][pos_i + old], en);
        pos_j+=K;
    }

    return en;
}

pub fn total_energy(system: &Vec<i8>, cplngs: &Vec<Vec<f32>>) -> f32 {
    let mut en:f32 = 0.0;
    let mut pos_i :usize = 0;
    for aa_i in system.iter() {
        let mut pos_j :usize = 0;
        for aa_j in system.iter() {
            en += cplngs[pos_i + *aa_i as usize][pos_j + *aa_j as usize];
            pos_j += K;
        }
        pos_i += K;
    }
    return en/2.0;
}

pub fn accumulate_counts(system: &Vec<i8>, counts: &mut Vec<Vec<f32>>) {
    let mut pos_i: usize = 0;
    for aa_i in system.iter() {
        let mut pos_j: usize = 0;
        for aa_j in system.iter() {
            counts[pos_i + *aa_i as usize][pos_j + *aa_j as usize] += 1.0;
            pos_j += K;
        }
        pos_i += K;
    }
}

pub fn print_matrix(m: &Vec<Vec<f32>>) {
    for row in m.iter() {
        for val in row.iter() { print!(" {:.3}",val); }
        println!("");
    }
}

pub fn isothermal_mc(system: &mut Vec<i8>, cplngs: &Vec<Vec<f32>>, inner_cycles: i32, outer_cycles: i32) -> Vec<Vec<f32>>{

    let mut counts: Vec<Vec<f32>> = vec![vec![0.0; N * K]; N * K];
    let mut rng = rand::thread_rng();
    let mut total_en:f32 = total_energy(&system,&cplngs);
    let mut n_obs: f32 = 0.0;
    for io in 0..outer_cycles {
        for ii in 0..inner_cycles {
            for is in 0..system.len() {
                let pos: usize = rng.gen_range(0..N);
                let new_aa: usize = rng.gen_range(0..K);
                let delta_en = delta_energy(system, cplngs, pos, system[pos] as usize, new_aa);
                if delta_en > 0.0 && (-delta_en).exp() < rng.gen_range(0.0..1.0) { continue; }
                system[pos] = new_aa as i8;
                total_en += delta_en;
            }
        }
        n_obs += 1.0;
        accumulate_counts(&system, &mut counts);
        // println!("{} {:?}",total_en,system);
    }
    counts.iter_mut().for_each(|el| el.iter_mut().for_each(|iel| *iel /= n_obs));
    return counts;
}

pub fn main() {
    let mut system: Vec<i8> = vec![0; N];
    let mut cplngs: Vec<Vec<f32>> = vec![vec![0.0; N * K]; N * K];

    init_1d(&mut cplngs);
    // println!("{:?}", &cplngs);

    total_energy(&system, &cplngs);
    // delta_energy(&system, &cplngs,0,0,1);
    let counts = isothermal_mc(&mut system, &cplngs,100,10000);
    // println!("{:?}",counts);
    print_matrix(&counts);
}