// Helper functions for supersingularity testing of Jacobians of genus-2 hyperelliptic curves over Fp via point multiplication (Algorithm 4), Cartier-Manin matrix, point counting and L-polynomial computation;
// Helper functions for superspeciality testing of Jaconians of genus-2 hyperelliptic curves over Fp2 via point multiplication (Algorithm 3) and Cartier-Manin matrix.
// Helper functions for conclusive sueprsingularity and superspeciality testing of Jacobians of genus-2 hyperelliptic curves over Fp^2 via pairings


// Comparison functions

function TestCartierSSing(p,J);
// TestCartierSSing(p::RngIntElt,J::JacHyp) -> BoolElt
// Returns true if the given Jacobian J of a genus-2 hyperelliptic curve over GF(p) has Cartier Representation C such that Cp * C = 0 where Cp is the matrix with entries of C raised to the p-power
// Cfr Lemma 1.1 in Ibukiyama-Katsura-Oort, 1986

    // assert IsPrime(p);
    // assert Characteristic(BaseField(J)) eq p;
    // assert Genus(Curve(J)) eq 2;

    C := Curve(J);
    CR := CartierRepresentation(C);
    CRp:= Parent(CR)!0;
    for ind in CartesianPower([1,2],2) do
        i := ind[1];
        j := ind[2];
        CRp[i][j] := CR[i][j]^p;
    end for;
    return CRp*CR eq Parent(CR)!0;
end function;

function TestPointCount(p,J,f);
// TestPointCount(p::RngIntElt,J::JacHyp,f::RngUPolElt) -> BoolElt
// Returns true if the given Jacobian J of a genus-2 hyperelliptic curve over GF(p) and its twist have the same number p^2 + b*p + 1 of points, for some b in {-2, -1, 0, 1, 2}, and otherwise false

    // assert IsPrime(p);
    // assert Characteristic(BaseField(J)) eq p;
    // assert Degree(f) in [5,6];
    // assert J eq Jacobian(HyperellipticCurve(f));
    
    ords := Orders(p);

    ord_J := Order(J);
    if not ord_J in ords then 
        return false;
    end if;

    Jt, _ := TwistedJacobian(p,f);
    ord_Jt := Order(Jt);
    if not ord_Jt in ords then 
        return false;
    end if;

    return ord_J eq ord_Jt;
end function;

function TestLPoly(p,J);
// TestLpoly(p::RngIntElt,J::JacHyp) -> BoolElt
// Returns true if the given Jacobian J of a genus-2 hyperelliptic curve over GF(p) has L-polynomial with coefficients [1, 0, b*p, 0, p^2] for some b

    // assert IsPrime(p);
    // assert Degree(f) in [5,6];
    // assert J eq Jacobian(HyperellipticCurve(f));

    C := Curve(J);
    Lpoly := LPolynomial(C);
    coeffs := Coefficients(Lpoly);

    // assert coeffs[1] eq 1;
    // coeffs[4] eq p*coeffs[2];

    if coeffs[2] ne 0 then 
        return false;
    elif coeffs[3] mod p ne 0 then
        return false;
    elif coeffs[5] ne p^2 then 
        return false;
    else
        return true;
    end if;
end function;

function TestCartierSSpecial(p,J)
// TestCartierSSpecial(p::RngIntElt,J::JacHyp) -> BoolElt
// Returns true if the given Jacobian J of a hyperelliptic curve in characteristic p has zero Cartier Representation C 
// Cfr Lemma 1.1 in Ibukiyama-Katsura-Oort, 1986

    // assert IsPrime(p);
    // assert Characteristic(BaseField(J)) eq p;

    C := Curve(J);
    CR := CartierRepresentation(C);
    if CR eq Parent(CR)!0 then 
        return true;
    else
        return false;
    end if;
end function;

function PairingVerification(p,J,f,n)
// PairingVerification(p::RngIntElt,J::JacHyp,f::RngUPolElt,n::RngIntElt) -> BoolElt
// Helper function for pairing-based supersingularity testing based on Section 5

    ords := Orders(p);

    if n eq ords[5] then // b = -2, n = (p+1)^2
        b := -2;
        m := GetExponent(p,n,b,f);
        // printf "b = %o\n", b;
        return PairingVerificationTypeMin2(J,m);
    elif n eq ords[4] then // b = -1, n = p^2 + p +1
        b := -1;
        m := GetExponent(p,n,b,f);
        // printf "b = %o\n", b;
        return PairingVerificationTypeMin1(J,m);
    elif n eq ords[3] then // b = 0, n = p^2 + 1
        b := 0;
        m := GetExponent(p,n,b,f);
        // printf "b = %o\n", b;
        return PairingVerificationType0(J,m);
    elif n eq ords[2] then // b = 1, n = p^2 - p + 1
        b := 1;
        m := GetExponent(p,n,b,f);
        // printf "b = %o\n", b;
        return PairingVerificationType1(J,m);
    elif n eq ords[1] then // b = 2, n = (p-1)^2
        b := 2;
        assert n eq (p-1)^2;
        m := GetExponent(p,n,b,f);
        // printf "b = %o\n", b;
        return PairingVerificationType2(J,m);
    else 
        return false;
    end if;
end function;

