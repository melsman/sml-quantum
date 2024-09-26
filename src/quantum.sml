
(* Gates *)

signature QUANTUM = sig
  type t                          (* Circuits *)
  val I            : t            (* Identity *)
  val X            : t            (* Pauli X *)
  val Y            : t            (* Pauli Y *)
  val Z            : t            (* Pauli Z *)
  val H            : t            (* Hadamard *)
  val C            : t -> t       (* Control *)
  val **           : t * t -> t   (* Tensor product *)
  val oo           : t * t -> t   (* Composition *)

  val ppC          : t -> string

  val check        : t -> int * int

  type ket
  val ket          : int list -> ket
  val ppK          : ket -> string

  type state
  val ppS          : state -> string
  val init         : ket -> state
  val eval         : t -> state -> state

  type dist = (ket*real) vector
  val ppDist       : dist -> string
  val measure      : state -> ket          (* Random measure *)
  val measurements : int -> state -> dist  (* Run many measurements *)

  type complex = Complex.complex
  type mat = complex Matrix.t
  val ppM          : mat -> string
  val sem          : t -> mat
end

structure Quantum :> QUANTUM = struct

infix |>
fun x |> f = f x

datatype t = I | X | Y | Z | H
           | Tensor of t * t
           | Compose of t * t
           | C of t

fun ppC t =
    case t of
        I => "I" | X => "X" | Y => "Y" | Z => "Z" | H => "H"
        | Tensor(t1,t2) => ppC t1 ^ " ** " ^ ppC t2
        | Compose(t1,t2) => ppC t1 ^ " oo " ^ ppC t2
                          | C t => "C(" ^ ppC t ^ ")"

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
type state = vec

fun one () = C.fromInt 1
fun zero () = C.fromInt 0

fun pow2 0 = 1
  | pow2 n = 2 * pow2(n-1)

type ket = int list
fun ket xs = xs
fun init (is: int list) : state =
    let val v = Vector.fromList is
        val i = Vector.foldli(fn (i,x,a) => x * pow2 i + a) 0 v
    in Vector.tabulate(pow2 (Vector.length v),
                       fn j => if i = j then one() else zero())
    end

fun toKet (n:int) (i:int) : ket =   (* state i \in [0;2^n-1] among total states 2^n, in binary *)
    let val s = StringCvt.padLeft #"0" n (Int.fmt StringCvt.BIN i)
    in CharVector.foldr (fn (#"1",a) => 1::a | (_,a) => 0::a) nil s
    end

fun matmul (t1:mat,t2:mat) : mat =
    M.matmul_gen C.* C.+ (zero()) t1 t2

fun tensor (t1: mat,t2:mat) : mat =
    M.tabulate(M.nRows t1 * M.nRows t2, M.nCols t1 * M.nCols t2,
               fn (r,c) =>
                  let val r1 = r div (M.nRows t2)
                      val c1 = c div (M.nCols t2)
                      val r2 = r mod (M.nRows t2)
                      val c2 = c mod (M.nCols t2)
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

fun ppS (v:vec) : string =
    let val v = Vector.map ppc v
        val sz = Vector.foldl (fn (e,a) => Int.max(a,size e)) 1 v
    in M.ppv sz (fn x => x) v
    end

fun ppK (v:ket) : string =
    "|" ^ implode (map (fn i => if i > 0 then #"1" else #"0") v) ^ ">"

fun log2 n = if n <= 1 then 0
             else 1 + log2 (n div 2)

fun for N f =
    let fun repeat i =
            if i >= N then ()
            else (f(); repeat (i+1))
    in repeat 0
    end

val gen = Random.newgen()

fun measure0 (s:state) : int =
    let fun square r = r * r
        val ps = Vector.map (square o Complex.mag) s
        val r : real = Random.random gen
    in Vector.foldl (fn (x,(i,a)) =>
                        let val a' = x+a
                        in if a' > r then (i,a')
                           else (i+1,a')
                        end) (0,0.0) ps
                    |> #1
    end

fun measure (s:state) : ket =
    toKet (log2 (Vector.length s))
          (measure0 s)

type dist = (ket * real) vector
fun ppDist (d:dist) : string =
    Vector.foldr (fn ((k,r),a) =>
                     (ppK k ^ " : " ^ Real.toString r) :: a) nil d
                 |> String.concatWith "\n"

fun measurements N s : dist =
    if N < 1 then raise Fail "measuments: expects positive N"
    else let val a = Array.array(Vector.length s, 0)
             val () = for N (fn () =>
                                let val i = measure0 s
                                in Array.update(a,i,Array.sub(a,i)+1)
                           end)
             val n = log2 (Vector.length s)
         in Vector.tabulate(Array.length a,
                            fn i => (toKet n i, real (Array.sub(a,i)) / real N))
         end

val $ = map (map Complex.fromInt)
val m1 = M.fromListList ($[[1,2,3],[4,5,6]])
val m2 = M.fromListList ($[[7,8],[9,10]])
val m3 = tensor (m1, m2)
val () = print (ppM m3 ^ "\n")

end


structure Test = struct

  open Quantum

  fun runSystem (c:t) k N =
      let val () = print ("Semantics of c = " ^ ppC c ^ ":\n")
          val () = print (ppM (sem c) ^ "\n")

          val s = init k
          val () = print ("Initial state s: " ^ ppS s ^ "\n")
          val s' = eval c s
          val () = print ("eval c " ^ ppK k ^ " = " ^ ppS s' ^ "\n")

          val () = print ("Distribution of " ^ Int.toString N ^ " measurements:\n")
          val () = print (ppDist(measurements N s') ^ "\n")
          val () = print "-----------------------------------\n"
      in ()
      end

  infix oo **

  val () = print "-----------------------------------\n"

  val () = runSystem H (ket[0]) 1000

  val () = runSystem (H oo H) (ket[0]) 1000

  val () = runSystem ((H oo H) ** I) (ket[0,1]) 1000

  val () = runSystem (I ** I) (ket[0,1]) 1000

end
