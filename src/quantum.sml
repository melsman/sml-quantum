
(* Gates *)

signature QUANTUM = sig
  type t
  val I : t             (* Identity *)
  val X : t             (* Pauli X *)
  val Y : t             (* Pauli Y *)
  val Z : t             (* Pauli Z *)
  val H : t             (* Hadamard *)
  val C : t -> t        (* Control *)
  val ** : t * t -> t   (* Tensor product *)
  val oo : t * t -> t   (* Composition *)

  val check : t -> int * int

  type complex = Complex.complex
  type mat = complex Matrix.t
  type vec = complex vector

  val ket : int list -> vec
  val eval : t -> vec -> vec
  val sem  : t -> mat

  val ppM : mat -> string
  val ppV : vec -> string
end

structure Quantum :> QUANTUM = struct

datatype t = I | X | Y | Z | H | Tensor of t * t | Compose of t * t | C of t

val op ** = Tensor
val op oo = Compose

fun check (t:t) : int * int =
    case t of
        Tensor (t1,t2) =>
        let val (i1,o1) = check t1
            val (i2,o2) = check t2
        in (i1*i2,o1*o2)
        end
      | Compose (t1,t2) =>
        let val (i1,o1) = check t1
            val (i2,o2) = check t2
        in if o1 <> i2 then raise Fail "Compose failure"
           else (i1,o2)
        end
      | C t' =>
        let val (i',o') = check t'
        in (1+i',1+o')
        end
      | _ => (1,1)

structure C = Complex
structure M = Matrix
type complex = C.complex
type mat = complex M.t
type vec = complex vector

fun one () = C.fromInt 1
fun zero () = C.fromInt 0

fun pow2 0 = 1
  | pow2 n = 2 * pow2(n-1)

fun ket (is: int list) : vec =
    let val v = Vector.fromList is
        val i = Vector.foldli(fn (i,x,a) => x * pow2 i + a) 0 v
    in Vector.tabulate(pow2 (Vector.length v),
                       fn j => if i = j then one() else zero())
    end

fun matmul (t1:mat,t2:mat) : mat =
    M.matmul_gen C.* C.+ (zero()) t1 t2

fun tensor (t1: mat,t2:mat) : mat =
    M.tabulate(M.nRows t1 * M.nRows t2, M.nCols t1 * M.nCols t2,
               fn (r,c) =>
                  let val r1 = r div (M.nRows t2)
                      val c1 = c div (M.nCols t2)
                      val r2 = r mod (M.nRows t2)
                      val c2 = c div (M.nCols t2)
                  in C.* (M.sub(t1,r1,c1), M.sub(t2,r2,c2))
                  end)

fun i () = C.fromIm 1.0

fun opI () = M.tabulate(2,2,fn(r,c) => if r = c then one() else zero())

fun opZ () = M.tabulate(2,2,fn(r,c) => if r <> c then zero()
                                       else if r = 0 then one()
                                       else C.~(one()))

fun opY () = M.tabulate(2,2,fn(r,c) => if r = c then zero()
                                       else if r = 0 then C.~(i())
                                       else i())

fun opX () = M.tabulate(2,2,fn(r,c) => if r = c then zero() else one())

fun opH () = M.tabulate(2,2,fn(r,c) => let val rsqrt2 = C.fromRe (1.0 / Math.sqrt 2.0)
                                       in if r = 1 andalso c = 1 then C.~ rsqrt2
                                          else rsqrt2
                                       end)

fun opC t = raise Fail "Control: not implemented"

fun sem (t:t) : mat =
    case t of
        I => opI()
      | X => opX()
      | Y => opY()
      | Z => opZ()
      | H => opH()
      | C t => opC t
      | Compose(t1,t2) => matmul(sem t1,sem t2)
      | Tensor(t1,t2) => tensor(sem t1,sem t2)

fun eval (t:t) (v:vec) : vec =
    M.matvecmul_gen C.* C.+ (zero()) (sem t) v

fun ppc (c:complex) =
    let val s = C.fmt (StringCvt.GEN(SOME 4)) c
    in if String.isSuffix ".0" s then
         String.extract(s,0,SOME(size s - 2))
       else s
    end

fun ppM (m:mat) : string =
    let val m = M.map ppc m
        val sz = foldl (fn (e,a) => Int.max(a,size e)) 1
                       (List.concat (M.listlist m))
    in M.pp sz (fn x => x) m
    end

fun ppV (v:vec) : string =
    let val v = Vector.map ppc v
        val sz = Vector.foldl (fn (e,a) => Int.max(a,size e)) 1 v
    in M.ppv sz (fn x => x) v
    end
end


structure Test = struct

  open Quantum

  val k = ket [0]

  infix oo **
  val system = H oo H

  val () = print ("sem :\n" ^ ppM (sem system) ^ "\n")

  val () = print ("eval S(" ^ ppV k ^ ") = " ^ ppV (eval system k) ^ "\n")
end
