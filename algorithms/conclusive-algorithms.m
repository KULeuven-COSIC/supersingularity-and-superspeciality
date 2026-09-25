// This file contains the pairing-based conclusive algorithms from section 5 of the associated paper.
// As in the paper, these differ per b-type of the Jacobian

function PairingVerificationType2(J,m)
// PairingVerificationType2(J::JacHyp,m::RngIntElt) -> BoolElt
// Pairing-based conclusive supersingularity testing for b = 2 based on Section 5, corollary 5.5
// Supersingular Jacobian(J) of genus-2 hyperelliptic curve over GF(p) with exponent m in [p-1, (p-1) div 2]

    p := Characteristic(BaseField(J));
    Fp := GF(p);
    Fp2 := GF(p^2);
    J2 := BaseExtend(J,Fp2);

    // Sample m-torsion points P1, P2 on J and Q1, Q2 on Jt until they give pairings of correct order
    // In this case, we first perform the test on J, and afterwards on J^t
    while true do
        P1 := Random(J); P2 := Random(J);
        if Order(WeilPairing(P1, P2, m)) eq (m div 2) then
            break;
        end if;
    end while;

    while true do
        R1 := Random(J2);
        R2 := Random(J2);
        Q1 := R1 - Pi(R1);
        Q2 := R2 - Pi(R2);

        if Order(WeilPairing(Q1, Q2, m)) eq (m div 2) then
            return true;
        else 
            continue;
        end if;
    end while;
end function;

function PairingVerificationTypeMin2(J,m)
// PairingVerificationTypeMin2(J::JacHyp,m::RngIntElt) -> BoolElt
// Pairing-based supersingularity testing for b=-2 based on Section 5, Corollary 5.5
// Supersingular Jacobian(J) of genus-2 hyperelliptic curve over GF(p) with exponent m in [p+1, (p+1) div 2]

    p := Characteristic(BaseField(J));
    Fp := GF(p);
    Fp2 := GF(p^2);
    J2 := BaseExtend(J,Fp2);

    while true do
        // Sample m-torsion points in J2
        list_P := [Random(J2) : i in [1..4]];

        // Compute pairings and discrete logarithms
        M, logM := WeilMatrix(p, list_P);
        
        // In practice, we perform this as Algorithm 5
        SF := SmithForm(logM);
        D := [ SF[i,i] : i in [1..4] | SF[i,i] ne 0 ];
        gcd := [ Integers() ! GCD(p+1, di) : di in D] ;
        fact_N := [ (p+1) div g : g in gcd ];
        N := &*fact_N;
        if N ge 8*p*(p^2 + 1) then 
            return true;
        else 
            continue;
        end if;
    end while;
end function;

function PairingVerificationType0(J,m)
// PairingVerificationType0(J::JacHyp,m::RngIntElt) -> BoolElt
// Pairing-based supersingularity testing for b=0 based on Section 5, Corollary 5.10
// Supersingular Jacobian(J) of genus-2 hyperelliptic curve over GF(p) with exponent m = p^2 +1

    p := Characteristic(BaseField(J));
    Fp := GF(p);
    Fp2 := GF(p^2);
    Fp4 := GF(p^4);

    J2 := BaseExtend(J,Fp2);
    J4 := BaseExtend(J2,Fp4);
    
    ok_ord := [m, m div 2];

    // for efficiency, we split the pairings between P1 and Q1 from P2 and Q2
    while true do
        // Sample m-torsion points P1 on J(Fp) 
        R1 := Random(J2);
        P1 := J4 ! (R1 + Pi(R1));

        // Sample m-torsion points Q1 on J(Fp2)
        R2 := Random(J4);
        Q1 := Pi(Pi(Pi(R2))) + p*Pi(Pi(R2)) - Pi(R2) - p*R2; // pi(Q1) = p*Q1

        // Compute pairings
        if Order(WeilPairing(P1,Q1,m)) in ok_ord then 
            break;
        end if;
    end while;

    while true do 
        // Sample m-torsion points P1 on J(Fp) and P2 on Jt(Fp)
        R1 := Random(J2);
        P2 := J4 ! (R1 - Pi(R1));

        // Sample m-torsion points Q1 on J(Fp2) and Q2 on Jt(Fp2)
        R2 := Random(J4);
        Q2 := Pi(Pi(Pi(R2))) - p*Pi(Pi(R2)) - Pi(R2) + p*R2; // pi(Q2) = -p*Q2

        // Compute pairings
        if Order(WeilPairing(P2,Q2,m)) in ok_ord then 
            return true;
        else 
            continue;
        end if;
    end while;
end function;


function PairingVerificationType1(J,m)
// PairingVerificationType1(p::RngIntElt,J::JacHyp,m::RngIntElt) -> BoolElt
// Pairing-based supersingularity testing for b=1 based on Section 5, Corollary 5.10
// Supersingular Jacobian(J) of genus-2 hyperelliptic curve over GF(p) with exponent m = p^2 -p + 1

    p := Characteristic(BaseField(J));
    Fp := GF(p);
    Fp2 := GF(p^2);
    Fp6 := GF(p^6);

    J2 := BaseExtend(J,Fp2);
    J6 := BaseExtend(J2,Fp6);

    if p mod 3 eq 1 then 
        ord_ok := m;
    else // p mod 3 eq 2 
        ord_ok := m div 3;
    end if;

    // For efficiency, we first compute the pairing between P1 and Q1, and then P2 and Q2
    while true do
        // Sample m-torsion point P1 on J(Fp)
        R1 := Random(J2);
        P1 := J6 ! (R1 + Pi(R1));

        // Sample m-torsion point Q1 on J(Fp6)
        R2 := (p+1)*Random(J6); // *(p+1) to make sure this is m-torsion
        Q1 := Pi(Pi(Pi(R2))) + p*Pi(Pi(R2)) - Pi(R2) - p*R2;

        if Order(WeilPairing(P1, Q1, m)) eq ord_ok then 
            break;
        end if;
    end while;

    while true do
        // Sample m-torsion points P2 on Jt(Fp)
        R1 := Random(J2);
        P2 := J6 ! (R1 - Pi(R1)) ;

        // Sample m-torsion points Q2 on J(Fp6)
        R2 := (p+1)*Random(J6); // *(p+1) to make sure this is m-torsion
        Q2 := Pi(Pi(Pi(R2))) - p*Pi(Pi(R2)) - Pi(R2) + p*R2;

        if Order(WeilPairing(P2, Q2, m)) eq ord_ok then 
            return true;
        else 
            continue;
        end if;
    end while;
end function;

function PairingVerificationTypeMin1(J,m)
// PairingVerificationType1(J::JacHyp,m::RngIntElt) -> BoolElt
// Pairing-based supersingularity testing for b=-1 based on Section 5, Corollary 5.10
// Supersingular Jacobian(J) of genus-2 hyperelliptic curve over GF(p) with exponent m = p^2 +p + 1

    p := Characteristic(BaseField(J));
    Fp := GF(p);
    Fp2 := GF(p^2);
    Fp6 := GF(p^6);

    J2 := BaseExtend(J,Fp2);
    J6 := BaseExtend(J2,Fp6);

    if p mod 3 eq 1 then 
        ord_ok := m div 3;
    else // p mod 3 eq 2
        ord_ok := m;
    end if;

    // For efficiency, we first compute the pairing between P1 and Q1, and then P2 and Q2
    while true do
        // Sample m-torsion point P1 on J(Fp)
        R1 := Random(J2);
        P1 := J6 ! (R1 + Pi(R1));

        // Sample m-torsion point Q1 on J(Fp3) 
        R2 := (p-1)*Random(J6); // *(p-1) to make sure this is m torsion 
        Q1 := Pi(Pi(Pi(R2))) + p*Pi(Pi(R2)) - Pi(R2) - p*R2;

        if Order(WeilPairing(P1, Q1, m)) eq ord_ok then 
            break;
        end if;
    end while;

    while true do
        // Sample m-torsion point P2 on Jt(Fp)
        R1 := Random(J2);
        P2 := J6 ! (R1 - Pi(R1));

        // Sample m-torsion point Q2 on J(Fp3) 
        R2 := (p-1)*Random(J6); // *(p-1) to make sure this is m torsion 
        Q2 := Pi(Pi(Pi(R2))) - p*Pi(Pi(R2)) - Pi(R2) + p*R2;

        if Order(WeilPairing(P2, Q2, m)) eq ord_ok then
            return true;
        else 
            continue;
        end if;
    end while;
end function;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//////   FOR CRYPTOGRAPHIC PRIMES p = c*2^f - 1 ////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

function PairingVerificationTypeMin2_crypto(J,m)
// PairingVerificationTypeMin2_crypto(J::JacHyp, m::RngIntElt) -> BoolElt
// Supersingularity testing based on Lemma 5.15 integrating 2-torsion shortcut
// Prime p = c*2^k -1 with (c odd and) c << 2^k

    p := Characteristic(BaseField(J));
    Fp := GF(p);

    k := Valuation(p+1, 2); // Largest k such that 2^k | (p+1)

    // Sample 2^k torsion points on J(Fp) 
    cof := (p+1) div (2^k);

    while true do
        P1 := cof * Random(J);
        P2 := cof * Random(J);

        // Order of P1
        k1 := 1;
        while true do 
            P1_double := 2*P1;
            if not P1_double eq J!0 then 
                k1 +:= 1;
                P1 := P1_double;
            else 
                break;
            end if;
        end while;

        // Order of P2
        k2 := 1;
        while true do 
            P2_double := 2*P2;
            if not P2_double eq J!0 then 
                k2 +:= 1;
                P2 := P2_double;
            else 
                break;
            end if;
        end while;
         
        bound := 3 + 3*Log(2,p);

        if k1 + k2 gt bound then 
            return true;
        else 
            continue;
        end if;
    end while;
end function;


function PairingVerificationCrypto(J,f,n)
// PairingVerificationCrypto(J::JacHyp,f::RngUPolElt,n::RngIntElt) -> BoolElt
// Pairing-based supersingularity testing improved for primes of the form p = c2^k -1 based on Lemma 5.13
// Only implemented for b = -2

    p := Characteristic(BaseField(J));
    b := -2;
    m := GetExponent(p,n,b,f);
    return PairingVerificationTypeMin2_crypto(J,m);
end function;