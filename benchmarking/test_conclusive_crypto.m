// Benchmarking conclusive supersingularity testing 
// for primes of the form p = c2^k -1
clear;

load "../helpers/arithmetic_helper_functions.m";
load "../algorithms/conclusive-algorithms.m";
load "../helpers/benchmark_helper_functions.m";

function CryptoPrimesBound(min_p, max_p)
    // CryptoPrimes(min_p::RngIntElt,max_p::RngIntElt) -> SeqEnum
    // Constructs a list of primes of the form p = c2^k -1 with c << 2^k, between min_p and max_p 
    P := [];
    for k in [9..800] do
        p := 5*2^k - 1;
        if p gt max_p then 
            break; 
        end if;
        q := (5*3)*2^k - 1;
        r := (7*5*3)*2^k - 1;
        s := (5*3^2)*2^k -1;
        L := [p,q,r,s];

        for x in L do    
            if IsPrime(x) and (x ge min_p) then 
                Append(~P,<x, Factorization(x+1)>); 
            end if; 
        end for;

    end for;
    return Sort(P);
end function;

samples_per_prime := 100; // Number of curves tested per prime
output_dat_file := "../benchmark_data/benchmark_data_conc_crypto.dat"; // Output file for TikZ / pgfplots

// Initialize / clear the data file with column headers
header := "p PairingCrypto_Avg Pairing_Avg \n";
Write(output_dat_file, header);

min_p := 29;
max_p := 2^(300);
CP := CryptoPrimesBound(min_p,max_p);
primes := [CP[i][1] : i in [1..#CP]];

printf "Starting benchmark across %o primes (%.2o bits to %.2o bits)\n", #primes, Log(2,Min(primes)), Log(2,Max(primes));
printf "Saving results to: %o\n\n", output_dat_file;

for p in primes do

    // Sample inputs
    input_Js := [];
    for i in [1..samples_per_prime] do 
        R<x> := PolynomialRing(GF(p));
        f := x^5 - 1;
        J := Jacobian(HyperellipticCurve(f));
        for j in [1..10] do 
            iso_Js := [ JJ : JJ in RichelotIsogenousSurfaces(J) | Type(JJ) eq Type(J) ];
            J := Random(iso_Js);
        end for;
        Append(~input_Js, <J, Polynomial(Curve(J))>);
    end for;

    printf "--------------------------------------------------\n";
    printf "Prime p = %o (%o curves available, testing %o)\n", p, #Set(input_Js), samples_per_prime;
    printf "--------------------------------------------------\n";

    // 1. Benchmark PairingVerificationCrypto
    t1 := Cputime();
    for j in [1..samples_per_prime] do
        item := input_Js[j];
        res := PairingVerificationCrypto(item[1],item[2],(p+1)^2);
    end for;
    t_tot_crypto := 10^6*Cputime(t1);
    avg_crypto := t_tot_crypto / samples_per_prime;

    // 2. Benchmark PairingVerification
    t2 := Cputime();
    for j in [1..samples_per_prime] do
        item := input_Js[j];
        res := PairingVerification(item[1],item[2],(p+1)^2); 
    end for;
    t_tot_pairing := 10^6*Cputime(t2);
    avg_pairing := t_tot_pairing / samples_per_prime;

    // Print summary to console
    printf "PairingCrypto: Avg = %o μs\n", avg_crypto;
    printf "Pairing: Avg = %o μs\n", avg_pairing;

    // Append formatted line to .dat file for plotting
    data_line := Sprintf("%o %o %o \n", p, avg_crypto, avg_pairing);
    PrintFile(output_dat_file, data_line);
end for;

printf "Benchmark complete. Results saved in '%o'\n", output_dat_file;
    
