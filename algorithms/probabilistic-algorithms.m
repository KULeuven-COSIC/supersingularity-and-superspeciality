// This file contains the probabilistic algorithms from section 4 of the associated paper.

function ProbabilisticSuperspecial(J)
// ProbabilisticSuperspecial(J::JacHyp,f::RngUPolElt) -> BoolElt
// Algorithm 3: probabilistic testing of superspeciality
    
    p := Characteristic(BaseField(J));
    C := Curve(J);
    twists := Twists(C);   

    for Crho in twists do
        Jrho := Jacobian(Crho);
        Jrho := BaseExtend(Jrho, GF(p^2));
        P := Random(Jrho);
        R := p*P;
        if R eq P or R eq -P then
            return true;
        end if; 
    end for;

    return false;
end function;

function IsSupersingular(J,f)
// IsSupersingular(J::JacHyp,f::RngUPolElt) -> BoolElt, RndIntElt
// Algorithm 4 (IsSupersingular)

    p := Characteristic(BaseField(J));

    repeat P := Random(J);
    until 2*P ne J!0 and 3*P ne J!0; // Avoid 2 and 3-torsion

    R := p*P;

    // (p+1)*P = 0?
    if R eq -P then  return true, (p+1)^2;  end if;
     
    Jt, _ := TwistedJacobian(p,f);
    repeat Q := Random(Jt);
    until 2*Q ne Jt!0 and 3*Q ne Jt!0; // Avoid 2 and 3-torsion

    // (p-1)*P = 0?
    if R eq P then 
        S := p*Q; 
        // (p-1)*Q = 0?
        if S eq Q then return true, (p-1)^2; end if;
    end if;

    T := p*R; // T = (p^2)*P
    if T eq -P then // (p^2 + 1)*P = 0
        n := p^2 + 1;
    elif T+P eq R then // (p^2 - p + 1)*P = 0
        n := p^2 - p + 1;
    elif T+P eq -R then // (p^2 + p + 1)*P = 0
        n := p^2 + p + 1;
    else 
        return false, 0;
    end if;

    // then verify on twist
    if n*Q eq Jt!0 then return true, n; end if;

    return false, 0;
end function;