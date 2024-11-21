# duplicated records are summed

    Code
      res$expr_mat
    Output
                                                      S1    S2
      RRJMTQGR_3_(N(F)(N(H(H(N))(H(N(H))))))       2e+07    NA
      NLEKJSTK_5_(N(N(H(H)(H(H)(H)))))             2e+07    NA
      NLEKJSTK_5_(N(N(H(H)(H))))                   2e+07    NA
      RRJMTQGR_3_(N(F)(N(H(H(N))(H(N(H(H)))))))       NA 2e+07
      RRJMTQGR_3_(N(F)(N(H(H(N(H)))(H(N(H(S)))))))    NA 2e+07
      JITQKR_1_(N(N(H(H)(H))))                        NA 2e+07

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

