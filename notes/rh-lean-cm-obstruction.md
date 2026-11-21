# Where Alpöge's Lean proof is tied to ℚ: the CM-container obstruction to ERH

*Notes from reading the source of the Lean project released at alpo.ge on 2026-11-20
(mirror: github.com/CocoSgt/riemann-hypothesis-lean, identical to the 14-part archive,
sha256 `9631b4f0…49c7fc`). I have read source text only; I have not built the project. Line
numbers are from that snapshot. Corrections welcome.*

## Summary

The note says the method "easily extends" to GRH for Dirichlet L-functions and that the Extended
Riemann Hypothesis is expected "with more time". The source shows exactly where the argument is
specific to characters of ℚ. That place is not a technicality.

1. **The character-specific interface is thin.** Only 1,137 of 175,187 files mention
   `DirichletCharacter`. The final step from primitive to all characters is a 30-line lemma.
2. **The per-character core runs over a container field `E` that must be CM.** `E` is built as a
   subfield of a cyclotomic field. It is proved CM *because it is abelian over ℚ*, and
   `[IsCMField E]` appears in 16,231 files.
3. **For a finite-order Hecke character of an imaginary quadratic field, of order ≥ 3, no such
   container can exist.** The field cut out by the character is dihedral over ℚ, so it is neither
   CM nor totally real. Every subfield of a CM field is CM or totally real.

So extending this architecture to ERH for imaginary quadratic fields, the case behind Bach's
class-group bound, class-group algorithms and so on, needs a new idea. Porting the code will not
be enough.

## 1. The interface

- `Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean:18`
  ```lean
  def ZeroFreeAll : Prop := ZeroFreeOn fun _ => True
  -- ZeroFreeOn Q := ∀ ψ : PChar, Q ψ → ∀ s, 1/2 < s.re → DirichletCharacter.LFunction ψ.χ s ≠ 0
  ```
  `PChar` is a primitive Dirichlet character.
- `Development/AutomorphicArithmetic/CharacterBoxes/Imprimitive.lean`: `ZeroFreeAll` implies the
  statement for all Dirichlet characters, by removing the finitely many Euler factors. This is
  elementary.
- `CharacterBoxes/EveryCharacter.lean:16` (`plugAll`): `ZeroFreeAll` follows from three inputs,
  `GLq.Green`, a `Hyps … DRS` family, and `VWin`. The character enters through `supplyChi ψ (hV ψ)`.
- `CharacterBoxes/ExistsPageN.lean:49`: `def VWin : Prop := ∀ ψ : PChar, VWinAt ψ 401`, a bound on
  Tate–Shafarevich lengths (`shaLen`) of Tate twists over number fields `E` with
  `IsChiField ψ e q d E`. This is the "weakened Vandiver" of the note's footnote 6.

## 2. The weld: `IsChiField` and CM-ness

`CharacterBoxes/ExistsPageN.lean:19`:
```lean
def IsChiField (ψ : PChar) (e q d : ℕ) (E : Type) [Field E] [NumberField E] : Prop :=
  ∃ (N : ℕ) (_ : NeZero N) (j : E →ₐ[ℚ] CyclotomicField N ℚ) (ζ : CyclotomicField N ℚ),
    2 ^ (e + 2) ∣ N ∧ ψ.N ∣ N ∧ q ∣ N ∧ d ∣ N ∧ IsPrimitiveRoot ζ N ∧
    ∀ τ a, τ ζ = ζ ^ a →
      ((∀ x : E, τ (j x) = j x) ↔
        (a ≡ 1 [MOD 2^(e+2)] ∨ a + 1 ≡ 0 [MOD 2^(e+2)]) ∧ ψ.val a = 1 ∧
          jacobiSym a q = 1 ∧ jacobiSym a d = 1)
```
So `E` is the subfield of ℚ(ζ_N) cut out by four conditions:
- the real 2-power cyclotomic layer,
- `ker ψ`,
- the quadratic character of q (q ≡ 3 mod 8),
- the quadratic character of d (d ≡ 7 mod 16).

That makes `E` the compositum of ℚ(ζ_{2^{e+2}})⁺, the field cut out by ψ, ℚ(√−q) and ℚ(√−d).
This matches the description in Manvir Jaswal's comment on Tao's blog.

`AutomorphicVarieties/CMFieldExponentialBounds/IsCMFieldOfIsChiField.lean:7`:
```lean
theorem UpperTheory.V.isCMField_of_isChiField :
    ∀ ψ e q d, q.Prime → q % 4 = 3 → ∀ E [Field E] [NumberField E], IsChiField ψ e q d E → IsCMField E
```
Its proof goes through `IsAbelianGalois ℚ E` (the embedding into a cyclotomic field) plus a Gauss
sum giving √−q ∈ E. Downstream, for example `AutomorphicVarieties/HeightAndVanishing/LayerLe.lean:29`,
`[IsCMField E]` is an instance argument. That is unsurprising: the theta lifts run between
*totally definite* unitary groups (note, step 8), which need a CM extension E/E⁺ with E⁺ totally
real.

## 3. Why this cannot be done for Hecke characters of imaginary quadratic fields

Let K be imaginary quadratic and ψ a character of Cl(K) of order n ≥ 3. Let L ⊂ H_K be the field
cut out by ψ.

- L is Galois over ℚ with Gal(L/ℚ) ≅ D_n. Complex conjugation acts on Gal(H_K/K) ≅ Cl(K) by
  inversion, and that action passes to the cyclic quotient of order n.
- If a Galois number field is CM, complex conjugation is central in its Galois group. In D_n with
  n ≥ 3, a reflection is not central. So L is not CM. It contains K, so it is not totally real either.
- Every subfield of a CM field is CM or totally real. **So no CM field contains L.**

Any container E for ψ that is analogous to `IsChiField` must contain L, the field through which
ψ factors, and must be CM. Both cannot hold. Over ℚ the problem never arises, because by
Kronecker–Weber every abelian field over ℚ is totally real or CM. That fact is exactly what
`isCMField_of_isChiField` uses.

The same problem affects finite-order Hecke characters of totally real fields F ≠ ℚ in general:
an abelian extension of F can have different complex conjugations at different real places, so
it need not lie in any CM field.

## What still extends with no new ideas

- Dedekind zeta functions of **abelian** number fields: they factor into Dirichlet L-functions.
- Quadratic (genus) characters of class groups: these are already Dirichlet characters.

## A possible second dependency on ℚ (not verified in code)

Jaswal reports that the minus part of the Tate–Shafarevich bound in `VWin` is controlled by
generalized Bernoulli numbers. That is the classical route for fields abelian over ℚ. For dihedral
fields, the analogous control would come from Rubin's main conjecture for imaginary quadratic
fields. There is a `GaloisLibrary/BernoulliNumbers` folder in the development, but I have not
traced this path.

## The GL₂ / Artin route is blocked too

Smallest test case: K = ℚ(√−23), h_K = 3, ψ a cubic class-group character. Its field L is the
Hilbert class field, the splitting field of x³ − x − 1, with Gal(L/ℚ) ≅ S₃. Then
L(s, ψ) = L(s, ρ) for the 2-dimensional representation ρ of S₃. This is the L-function of the
weight-one form η(z)η(23z). Also ζ_F = ζ · L(s, ρ) for the cubic field F = ℚ[x]/(x³ − x − 1),
so the target is equivalent to RH for ζ_F.

1. **L(s, ρ) is not a factor of ζ_E for any CM field E.** ρ occurs in Ind_E^ℚ 1 only if the image
   of Gal(Ẽ/E) in S₃ fixes a nonzero vector of ρ. That forces the image to be trivial or of
   order 2, so E ⊇ L or E ⊇ a conjugate of F. F has signature (1,1), so neither F nor L is CM or
   totally real. Neither can lie in a CM field.
2. **Changing the base field to K (which is CM) does not help.** The analogous container would be
   E/K abelian with ψ among its characters, so E ⊇ L again.
3. **The root cause is definiteness, not unitary groups as such.** The theta lifts run between
   totally definite groups. For a semisimple group, being compact at an infinite place is only
   possible at a real place, since at a complex place its points form a complex Lie group, which
   is noncompact. So any definite version needs a totally real field of definition, which for
   unitary groups means a CM container.

So extending this architecture beyond Dirichlet characters appears to need an indefinite version
of the "theta lifts mimic multiplicity" mechanism, on non-compact groups and Shimura varieties.
That would be new mathematics, not a port.
