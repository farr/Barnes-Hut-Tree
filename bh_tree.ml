(** A Barnes-Hut tree functor.  This code assumes that you are
    constructing a tree in three dimensions. *)

(*  bh_tree.ml: Barnes-Hut tree functor.
    Copyright (C) 2009 Will M. Farr <w-farr@northwestern.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*)

(** The type signature for the input argument to [Make]. *)
module type BODY = sig
  type b
        (** Body type. *)

  val q : b -> float array
      (** Get the position of the body. *)

  val summary_body : b array -> b
      (** [summary_body bs] constructs a body that summarizes the
          properties of the given [bs].  For example, it may construct
          a body representing the monopole moment (center of mass) of
          the given bodies.  Each cell in the tree contains a summary
          body that can be used to approximate the gravitational
          interaction with the bodies contained in that cell. *)

  val null_body : b
      (** A "null" value for bodies, such that it does not change the
          outcome of summary_body when included in the body array for
          that function. *)

end

(** Barnes-Hut trees. *)
module type BH_TREE = sig
  type b
        (** Bodies. *)

  type t = 
    | Empty (** Empty tree. *)
    | Body of b (** Singleton body. *)
    | Cell of float array * float array * b * t array (** Cell with bounds, summary body, and sub-trees. *)
          (** Trees. *)

  val empty : t 
      (** The empty tree. *)

  val in_bounds : float array -> float array -> float array -> bool
      (** [in_bounds q lower upper] checks whether [lower.(i) <= q.(i)
          < upper.(i)] for all elements of the arrays.  Note the
          inclusive lower- and exclusive upper-bound. *)

  val tree_of_body_list : b list -> t
      (** Equivalent to [List.fold_left (fun t b -> insert b t) Empty
          bs], but more efficient. *)

  val tree_of_body_array : b array -> t
      (** Equivalent to [Array.fold_right insert bs Empty], but more
          efficient. *)

  val insert : b -> t -> t
      (** Insert a body into the tree. *)

  val remove : b -> t -> t 
      (** Remove a body from the tree.  [remove b t] will only remove
          [b] from [t] if a body that is [==] to [b] is stored in the
          tree (i.e. [remove] does not compare bodies based on
          position). *)

  val bodies : t -> b list
      (** Returns a list of all bodies in the tree. *)

  val fold : (b -> 'a -> 'a) -> t -> 'a -> 'a
      (** [fold kons t knil] applies [kons] to each body in [t] (in
          unspecified order), using the result of each application as
          the second argument for the next.  [knil] is used as the
          second argument for the first application.  

          For example, the [bodies] function could be written [let
          bodies t = fold (fun b l -> b :: l) t \[\]] *)

  val fold_approx : (t -> bool) -> (b -> 'a -> 'a) -> t -> 'a -> 'a
      (** [fold_approx as_body kons t knil] is like [fold] except that
          the [as_body] is applied to all cells in [t].  If [as_body]
          returns [true], the summary body for that cell is used as
          the first argument to [kons], and the walk of that branch of
          the tree stops.  The function [fold_approx] is useful for
          algorithms that want to use the summary body as an
          approximation for all the bodies contained within a cell
          when some condition is met. *)

end

(** Functor to produce a Bh_tree module given a Body module. *)
module Make(B : BODY) : BH_TREE with type b = B.b = struct
  type b = B.b

  type t = 
    | Empty
    | Body of b
    | Cell of float array * float array * b * t array

  let empty = Empty

  let in_bounds q lower upper = 
    let inb = ref true and 
        i = ref 0 in 
    while !inb && (!i < 3) do 
      let qi = q.(!i) in 
      if qi < lower.(!i) || qi >= upper.(!i) then 
        inb := false;
      incr i
    done;
    !inb

(*   let split_bounds lower upper =  *)
(*     let split_bds = Array.make 8 (lower,upper) in  *)
(*     let mid = Array.make 3 0.0 in  *)
(*     for i = 0 to 2 do *)
(*       mid.(i) <- 0.5*.(lower.(i) +. upper.(i)) *)
(*     done; *)
(*     for i = 0 to 7 do  *)
(*       let ilower = Array.make 3 0.0 and  *)
(*           iupper = Array.make 3 0.0 in  *)
(*       for j = 0 to 2 do  *)
(*         if i land (1 lsl j) > 0 then begin *)
(*           (\* high in this dimension *\) *)
(*           ilower.(j) <- mid.(j); *)
(*           iupper.(j) <- upper.(j) *)
(*         end else begin *)
(*           ilower.(j) <- lower.(j); *)
(*           iupper.(j) <- upper.(j) *)
(*         end *)
(*       done; *)
(*       split_bds.(i) <- (ilower,iupper) *)
(*     done; *)
(*     split_bds *)

  let incorporate_body b lower upper = 
    let q = B.q b in 
    if in_bounds q lower upper then 
      ()
    else begin
      for i = 0 to 2 do 
        if q.(i) < lower.(i) then lower.(i) <- q.(i);
        if q.(i) >= upper.(i) then upper.(i) <- q.(i) +. 100.0*.epsilon_float*.(abs_float q.(i))
      done
    end

  let expand_bounds lower upper = 
    for i = 0 to 2 do 
      let delta = upper.(i) -. lower.(i) in 
      lower.(i) <- lower.(i) -. 0.5*.delta;
      upper.(i) <- upper.(i) +. 0.5*.delta
    done
    
  let sub_bounds_index q lower upper = 
    let i = ref 0 in 
    for j = 0 to 2 do 
      let lj = lower.(j) and 
          uj = upper.(j) and 
          qj = q.(j) in 
      let mj = 0.5*.(lj +. uj) in 
      if qj >= mj then i := !i + (1 lsl j)
    done;
    !i        

  let sub_bounds i lower upper = 
    let slow = Array.make 3 0.0 and 
        shigh = Array.make 3 0.0 in 
    for j = 0 to 2 do 
      if i land (1 lsl j) > 0 then begin
        slow.(j) <- 0.5*.(lower.(j) +. upper.(j));
        shigh.(j) <- upper.(j)
      end else begin
        slow.(j) <- lower.(j);
        shigh.(j) <- 0.5*.(lower.(j) +. upper.(j))
      end
    done;
    (slow, shigh)

  let bounds_of_body_list bs = 
    let lower = Array.make 3 infinity and 
        upper = Array.make 3 neg_infinity in 
    List.iter (fun b -> incorporate_body b lower upper) bs;
    expand_bounds lower upper;
    (lower,upper)

  let bounds_of_body_array bs = 
    let lower = Array.make 3 infinity and 
        upper = Array.make 3 neg_infinity and 
        n = Array.length bs in 
    for i = 0 to n - 1 do 
      incorporate_body bs.(i) lower upper
    done;
    expand_bounds lower upper;
    (lower,upper)

  let body_of_tree = function 
    | Empty -> B.null_body
    | Body(b) -> b
    | Cell(_,_,b,_) -> b

  let rec insert_w_bounds b lower upper = function 
    | Empty -> Body(b)
    | Body(b2) -> 
        let t = Cell(lower, upper, B.summary_body [|b; b2|], Array.make 8 Empty) in 
        insert_w_bounds b lower upper
          (insert_w_bounds b2 lower upper t)
    | Cell(lower, upper, _, sub_ts) -> 
        let q = B.q b in 
        let i = sub_bounds_index q lower upper in 
        let (sub_lower, sub_upper) = sub_bounds i lower upper in 
        let sub_ts = Array.copy sub_ts in 
        sub_ts.(i) <- insert_w_bounds b sub_lower sub_upper sub_ts.(i);
        let sum_b = B.summary_body (Array.map body_of_tree sub_ts) in 
        Cell(lower,upper,sum_b,sub_ts)

  let rec bodies = function 
    | Empty -> []
    | Body(b) -> [b]
    | Cell(_,_,_,sub_ts) -> 
        let bs = ref [] in 
        for i = 0 to 7 do 
          bs := (bodies sub_ts.(i)) @ !bs
        done;
        !bs

  let tree_of_body_array bs = 
    let (lower, upper) = bounds_of_body_array bs in 
    let t = ref Empty in 
    for i = 0 to Array.length bs - 1 do 
      t := insert_w_bounds bs.(i) lower upper !t
    done;
    !t

  let tree_of_body_list bs = 
    let (lower, upper) = bounds_of_body_list bs in 
    List.fold_left
      (fun t b -> insert_w_bounds b lower upper t)
      Empty
      bs

  let insert b = function 
    | Empty -> Body(b)
    | Body(b2) -> tree_of_body_array [|b; b2|]
    | Cell(lower, upper, _, _) as t -> 
        if in_bounds (B.q b) lower upper then 
          insert_w_bounds b lower upper t
        else
          let bs = b :: (bodies t) in 
          tree_of_body_list bs

  let all_empty ts = 
    let all_empty = ref true in 
    for i = 0 to 7 do 
      match ts.(i) with 
      | Empty -> ()
      | _ -> all_empty := false
    done;
    !all_empty

  let one_body ts = 
    let nbody = ref 0 and 
        any_cell = ref false in 
    for i = 0 to 7 do 
      match ts.(i) with 
      | Cell(_,_,_,_) -> any_cell := true
      | Body(_) -> incr nbody
      | _ -> ()
    done;
    not (!any_cell) && !nbody = 1

  let extract_body ts = 
    let b = ref B.null_body in 
    for i = 0 to 7 do 
      match ts.(i) with 
      | Body(bb) -> b := bb
      | _ -> ()
    done;
    !b

  let rec remove b = function 
    | Empty -> raise (Failure "remove: body not found in tree")
    | Body(b2) as t -> if b == b2 then Empty else t
    | Cell(lower, upper, _, sub_ts) -> 
        let q = B.q b in 
        if in_bounds q lower upper then 
          let i = sub_bounds_index q lower upper in 
          let new_sub = 
            Array.mapi (fun j sub -> if j = i then remove b sub else sub) sub_ts in 
          if all_empty new_sub then
            Empty
          else if one_body new_sub then 
            Body(extract_body new_sub)
          else
            let sb = B.summary_body (Array.map body_of_tree new_sub) in 
            Cell(lower, upper, sb, new_sub)
        else
          raise (Failure "remove: body not found in tree")

  let rec fold kons t knil = 
    match t with 
    | Empty -> knil
    | Body(b) -> kons b knil
    | Cell(_,_,_,sub_ts) -> 
        let res = ref knil in 
        for i = 0 to 7 do 
          res := fold kons sub_ts.(i) !res
        done;
        !res

  let rec fold_approx as_body kons t knil = 
    match t with 
    | Empty -> knil
    | Body(b) -> kons b knil
    | Cell(_,_,sb,sub_ts) as c -> 
        if as_body c then 
          kons sb knil
        else begin
          let res = ref knil in 
          for i = 0 to 7 do 
            res := fold_approx as_body kons sub_ts.(i) !res
          done;
          !res
        end
          

end
