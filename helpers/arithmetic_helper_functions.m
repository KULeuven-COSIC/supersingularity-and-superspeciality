// Helper functions for working with Jacobians over Fp

function Pi(P)
    p := Characteristic(BaseField(Parent(P)));
    Fp := GF(p);   
    return Frobenius(P, Fp);
end function;

function RandomHypCurvePoly(p)
// RandomHypCurvePoly(p::RngIntElt) -> RngUPolElt
// Constructs a random polynomial of a a genus-2 hyperellipitc curve over GF(p)

    // assert IsPrime(p);

    Fp := GF(p);
    R<x> := PolynomialRing(Fp);

    // Coin flip: degree 5 or 6
    deg := 5 + Random(1);

    // Random squarefree monic polynomial of the chosen degree
    f := R!0;
    repeat
        coeffs := [Random(Fp) : i in [1..deg]];
        f := x^deg + &+[coeffs[i]*x^(i-1) : i in [1..deg]];
    until IsSquarefree(f);
    return f;
end function;

function RandomJacobian(p)
// RandomJacobian(p::RngIntElt) -> JacHyp,RngUPolElt
// Constructs a random Jacobian of a genus-2 hyperellipitc curve over GF(p)

    f := RandomHypCurvePoly(p);
    C  := HyperellipticCurve(f);
    J  := Jacobian(C);

    return J, f;
end function;

function TwistedJacobian(p,f)
// TwistedJacobian(p::RngIntElt,f::RngUPolElt) -> JacHyp
// Constructs a quadratic twist of the Jacobian of the genus-2 hyperellipitc curve over GF(p) defined by y^2 = f

    // assert IsPrime(p);
    // assert Degree(f) in [5,6];

    Fp := GF(p);
    // Quadratic twist: find a non-square d in Fp, then y^2 = d*f(x)
    repeat
        d := Random(Fp);
    until d ne 0 and not IsSquare(d);

    Ct := HyperellipticCurve(d * f);
    return Jacobian(Ct), d;
end function;

function Orders(p)
// Orders(p::RngIntElt) -> SeqEnum
// Constructs the list of possible orders of supersingular Jacobians over GF(p), i.e. [p^2 + b*p + 1] for b in [-2,-1,0,1,2]

    p2 := p*p;
    p2p1 := p2 + 1;
    twop := p+p;
    ords := [p2p1 - twop, p2p1 - p, p2p1, p2p1 + p, p2p1 + twop];
    return ords;
end function;

function GetExponent(p,n,b,f);
// GetExponent(p::RngIntElt,b::RngIntElt,f::RngUPolElt) -> RngIntElt
// Returns the exponent of a supersingular Jacobian over GF(p) of a genus-2 hyperelliptic curve defined by f, with number of points n = p^2 -b*p + 1

    if b in [-1,0,1] then 
        return n;
    end if;

    if b eq -2 then 
        if p mod 4 eq 1 then 
            return p+1;
        else // p mod 4 eq 3 
            deg_f := Degree(f);
            fact_f := Factorization(f); 
            len_fact := #fact_f; // f separable, all factors distinct
            if len_fact lt deg_f -1 then 
                return p+1;
            else 
                return (p+1) div 2;
            end if;
        end if;
    end if;

    if b eq 2 then 
        if p mod 4 eq 1 then 
            deg_f := Degree(f);
            fact_f := Factorization(f); 
            len_fact := #fact_f; // f separable, all factors distinct
            if len_fact lt deg_f -1 then 
                return p-1;
            else 
                return (p-1) div 2;
            end if;
        else // p mod 4 eq 3
            return p-1;
        end if;
    end if;
end function;

function WeilMatrix(p,P)
// WeilMatrix(p::IntRngElt,P::SeqEnum) -> AlgMatElt, AlgMatElt
// Given a list of four points of order p+1 P, returns the matrix of Weil pairings of points in the list and the matrix of dicrete logarithms of such pairings with base a primitive (p+1)th root of unity
    
    Fp2 := BaseRing(Parent(P[1]));

    // Primitive (p+1)-th root of unity 
    zeta := PrimitiveElement(Fp2)^(p-1);

    M := ZeroMatrix(Fp2,4,4);
    for i in [1..4] do
        M[i][i] := 1;
        for j in [i..4] do
            M[i][j] := WeilPairing(P[i], P[j], p+1);
            M[j][i] := 1/M[i][j];
	end for;
    end for;

    Zm := ResidueClassRing(p+1);
    logM := ZeroMatrix(Zm, 4, 4);
    for i in [1..4] do
        logM[i][i] := 0;
        for j in [i..4] do
            logM[i][j] := Zm!Log(zeta, M[i][j]); 
	    logM[j][i] := -logM[i][j];
        end for;
    end for;

    return M, logM;
end function;