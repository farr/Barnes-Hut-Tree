open OUnit

let tests = "all bh_tree tests" >:::
  ["bh_tree.ml tests" >: Bh_tree_test.tests]

let _ = Random.self_init ()

let _ = 
  let results = run_test_tt_main tests in 
  let nfail = 
    List.fold_left
      (fun nfail test -> 
        match test with 
        | RSuccess(_) -> nfail
        | _ -> nfail + 1)
      0
      results in 
  exit nfail
