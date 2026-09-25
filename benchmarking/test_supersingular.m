// Benchmarking probabilistic supersingularity testing 
clear;

load "../helpers/arithmetic_helper_functions.m";
load "../algorithms/probabilistic-algorithms.m";
load "../algorithms/conclusive-algorithms.m";
load "../helpers/benchmark_helper_functions.m";
load "../helpers/ss_curves_29_211.m";

samples_per_prime := 1000; // Number of curves tested per prime
output_dat_file := "../benchmark_data/benchmark_data_ssing_prob.dat"; // Output file for benchmarks and TikZ / pgfplots
output_dat_file_falseneg := "../benchmark_data/false_neg_ssing_prob.dat"; // Output file for false negatives rates

// Initialize / clear the data file with column headers
header := "p PointMult_Avg PointCount_Avg\n";
Write(output_dat_file, header);

// header_falseneg := "p FalseNeg FalseNeg%%\n";
// Write(output_dat_file_falseneg, header_falseneg);

primes := Sort([p : p in Keys(SS_CURVES_BY_PRIME)]);

printf "Starting benchmark across %o primes (%.2o bits to %.2o bits)\n", #primes, Log(2,Min(primes)), Log(2,Max(primes));
printf "Saving results to: %o\n\n", output_dat_file;

for p in primes do
    curves_data := SS_CURVES_BY_PRIME[p];
    // samples_per_prime := Min(samples_per_prime, #curves_data);

    printf "--------------------------------------------------\n";
    printf "Prime p = %o (%o curves available, testing %o)\n", p, #curves_data, samples_per_prime;
    printf "--------------------------------------------------\n";

    // Define ring context for 'eval'
    Fp<x> := PolynomialRing(GF(p));

    // Sample inputs
    input_polys := [];
    for i in [1..samples_per_prime] do
        entry := Random(curves_data);
        f := eval entry[1];
        Append(~input_polys, f);
    end for;

    // 1. Benchmark IsSupersingular
    input_Js_polys := [ < Jacobian(HyperellipticCurve(f)), f > : f in input_polys ];
    results_point_mult := [];
    t1 := Cputime();
    for item in input_Js_polys do
        // J := item[1]; f := item[2];
        res, _ := IsSupersingular(item[1], item[2]);
        // Append(~results_point_mult,res);
    end for;
    t_tot_point_mult := 10^6*Cputime(t1);
    avg_point_mult := t_tot_point_mult / samples_per_prime; 

    // // 2. Benchmark TestLPoly
    // input_Js_polys := [ < Jacobian(HyperellipticCurve(f)), f > : f in input_polys ]; // Re-initialising Jacobians ensures no cached data is used
    // t2 := Cputime();
    // for item in input_Js_polys do
    //     // J := item[1]; f := item[2];
    //     _ := TestLPoly(p, item[1]);
    // end for;
    // t_tot_Lpoly := 10^6*Cputime(t2);
    // avg_Lpoly := t_tot_Lpoly / samples_per_prime;

    // 3. Benchmark TestPointCount
    input_Js_polys := [ < Jacobian(HyperellipticCurve(f)), f > : f in input_polys ]; // Re-initialising Jacobians ensures no cached data is used
    t3 := Cputime();
    for item in input_Js_polys do
        // J := item[1]; f := item[2];
        _ := TestPointCount(p, item[1], item[2]);
    end for;
    t_tot_point_count := 10^6*Cputime(t3);
    avg_point_count := t_tot_point_count / samples_per_prime;

    // // 4. Benchmark TestCartier
    // input_Js_polys := [ < Jacobian(HyperellipticCurve(f)), f > : f in input_polys ]; // Re-initialising Jacobians ensures no cached data is used
    // results_cartier := [];
    // t4 := Cputime();
    // for item in input_Js_polys do
    //     // J := item[1]; f := item[2];
    //     res := TestCartierSSing(p, item[1]);
    //     // Append(~results_cartier, res);
    // end for;
    // t_tot_cartier := 10^6*Cputime(t4);
    // avg_cartier := t_tot_cartier / samples_per_prime;

    // Print summary to console
    printf "PointMult : Avg = %o μs\n", avg_point_mult;
    // printf "LPoly     : Avg = %o μs\n", avg_Lpoly;
    printf "PointCount: Avg = %o μs\n", avg_point_count;
    // printf "Cartier   : Avg = %o μs\n\n", avg_cartier;

    // Append formatted line to .dat file for plotting
    data_line := Sprintf("%o %o %o \n", p, avg_point_mult, avg_point_count);
    PrintFile(output_dat_file, data_line);

    // // Checking no false negatives
    // false_neg_p := 0;
    // for j in [1..samples_per_prime] do 
    //     assert results_cartier[j];
    //     if results_point_mult[j] ne results_cartier[j] then 
    //         false_neg_p +:= 1;
    //         print "false_neg : ", input_Js_polys[j];
    //     end if;
    // end for;
    // printf "False negatives : %o\n", false_neg_p;
    // data_line_falseneg := Sprintf("%o %o %o\n", p, false_neg_p, false_neg_p * 100.000 / n);
    // PrintFile(output_dat_file_falseneg, data_line_falseneg);

end for;

printf "Benchmark complete. Results saved in '%o'\n", output_dat_file;