# duplicated records are summed

    Code
      res$expr_mat
    Output
                             S1    S2
      RRJMTQGR_3_H4N4F1   2e+07    NA
      NLEKJSTK_5_H5N2     2e+07    NA
      NLEKJSTK_5_H3N2     2e+07    NA
      RRJMTQGR_3_H5N4F1      NA 2e+07
      RRJMTQGR_3_H5N4F1A1    NA 2e+07
      JITQKR_1_H3N2          NA 2e+07

# differ_a_g setting to FALSE changes A to S

    Code
      res$glycan_graphs[[5]]
    Output
      Glycan Graph (NE)
      F: 1, H: 5, N: 4, S: 1
      ------------------
      N
      ├─F
      └─N
        └─H
          ├─H
          │ └─N
          │   └─H
          └─H
            └─N
              └─H
                └─S

# read_pglyco3 uses 'corrected' columns

    Code
      res_mono$expr_mat
    Output
                             S1    S2
      RRJMTQGR_3_H4N4A1   2e+07    NA
      NLEKJSTK_5_H3N4A1   2e+07    NA
      NLEKJSTK_5_H5N2A1   2e+07    NA
      RRJMTQGR_3_H5N4F1      NA 1e+07
      RRJMTQGR_3_H5N4F1A1    NA 1e+07
      JITQKR_1_H3N2          NA 1e+07

