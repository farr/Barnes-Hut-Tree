open OUnit

module B = struct 
  type b = {m : float;
            q : float array;
            p : float array}

  let make_random () = 
    {m = Random.float 1.0;
     q = Array.init 3 (fun _ -> Random.float 1.0 -. 0.5);
     p = Array.init 3 (fun _ -> Random.float 1.0 -. 0.5)}

  let q b = b.q

  let null_body = {m = 0.0; q = Array.make 3 0.0; p = Array.make 3 0.0}

  let summary_body bs = 
    let mtot = ref 0.0 and 
        qcom = Array.make 3 0.0 and 
        ptot = Array.make 3 0.0 in 
    for i = 0 to Array.length bs - 1 do  
      let {m = m; q = q; p = p} = bs.(i) in 
      mtot := !mtot +. m;
      for i = 0 to 2 do 
        qcom.(i) <- qcom.(i) +. m*.q.(i);
        ptot.(i) <- ptot.(i) +. p.(i)
      done
    done;
    for i = 0 to 2 do 
      qcom.(i) <- qcom.(i) /. !mtot
    done;
    {m = !mtot;
     q = qcom;
     p = ptot}
end
open B

module Bht = Bh_tree.Make(B)
open Bht

let afor_all pred arr = 
  let n = Array.length arr in 
  let i = ref 0 and 
      all = ref true in 
  while !i < n && !all do 
    if not (pred arr.(!i)) then all := false;
    incr i
  done;
  !all

let rec test_tree_valid = function 
  | Empty -> true
  | Body(_) -> true
  | Cell(lower, upper, _, sub_ts) as t -> 
      let bs = bodies t in 
      let contained = List.for_all (fun b -> in_bounds b.q lower upper) bs in 
      contained && (afor_all test_tree_valid sub_ts)

let test_array_list () = 
  let bs = Array.init 1000 (fun _ -> make_random ()) in 
  let bsl = Array.to_list bs in 
  let ta = tree_of_body_array bs and 
      tl = tree_of_body_list bsl in 
  assert_bool "array and list trees differ" (ta = tl)

let test_valid () = 
  let bs = Array.init 10000 (fun _ -> make_random ()) in 
  let t = tree_of_body_array bs in 
  assert_bool "some bodies not contained in their cell" (test_tree_valid t);
  assert_bool "not enough bodies in tree" (List.length (bodies t) = 10000)

let test_remove () = 
  let bs = Array.init 1000 (fun _ -> make_random ()) in 
  let b = bs.(Random.int 1000) in 
  let t = tree_of_body_array bs in 
  let tr = remove b t in 
  assert_bool "didn't remove a body" (List.length (bodies tr) = 999);
  assert_bool "didn't remove correct body" 
    (List.for_all (fun tb -> not (tb == b)) (bodies tr));
  let t' = insert b tr in 
  assert_bool "remove then insert not identity" (t = t')

let test_remove_to_empty () = 
  let bs = Array.init 1000 (fun _ -> make_random ()) in 
  let t = tree_of_body_array bs in 
  let loop_t = ref t in 
  for i = 0 to 999 do 
    loop_t := remove bs.(i) !loop_t
  done;
  assert_bool "didn't remove all bodies and find Empty" 
    (match !loop_t with 
    | Empty -> true
    | _ -> false)

let test_potential () = 
  let bs = Array.init 1000 (fun _ -> {(make_random ()) with m = 1e-3}) in 
  let distance q1 q2 = 
    let d2 = ref 0.0 in 
    for i = 0 to 2 do 
      let dx = q1.(i) -. q2.(i) in 
      d2 := !d2 +. dx*.dx
    done;
    !d2 in 
  let v = ref 0.0 in 
  for i = 0 to 999 do 
    let b = bs.(i) in 
    for j = i+1 to 999 do 
      let b2 = bs.(j) in 
      v := !v -. b.m*.b2.m/.(distance b.q b2.q)
    done
  done;
  let v_exact = !v in 
  let t = ref (tree_of_body_array bs) and 
      v = ref 0.0 in 
  for i = 0 to 999 do 
    let b = bs.(i) in 
    t := remove b !t;
    v := !v +. 
        (fold_approx 
           (function 
             | Cell(lower,upper,sb,_) -> 
                 let r = distance lower upper and 
                     d = distance b.q sb.q in 
                 r/.d < 0.1
             | _ -> false)
           (fun tb vsum -> 
             vsum -. b.m*.tb.m/.(distance b.q tb.q))
           !t
           0.0)
  done;
  assert_equal ~cmp:(cmp_float ~epsabs:0.0 ~epsrel:1e-2) 
    ~printer:string_of_float v_exact !v
                   
let tests = "bh_tree.ml tests" >:::
  ["array and list constructors" >:: test_array_list;
   "basic validity tests" >:: test_valid;
   "test remove" >:: test_remove;
   "test remove to empty" >:: test_remove_to_empty;
   "potential calculation" >:: test_potential]
