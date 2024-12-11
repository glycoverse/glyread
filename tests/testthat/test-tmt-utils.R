test_that("mz_to_tmt_channel works for TMT-6plex", {
  mz <- c(126.127726, 127.124761, 128.134436, 129.131471, 130.141145, 131.138180)
  expected <- c("126", "127", "128", "129", "130", "131")
  expect_equal(mz_to_tmt_channel(mz, "TMT-6plex"), expected)
})


test_that("mz_to_tmt_channel works for TMT-10plex", {
  mz <- c(126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
          129.131471, 129.137790, 130.134825, 130.141145, 131.138180)
  expected <- c("126", "127N", "127C", "128N", "128C",
                "129N", "129C", "130N", "130C", "131")
  expect_equal(mz_to_tmt_channel(mz, "TMT-10plex"), expected)
})


test_that("mz_to_tmt_channel works for TMT-11plex", {
  mz <- c(126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
          129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144499)
  expected <- c("126", "127N", "127C", "128N", "128C",
                "129N", "129C", "130N", "130C", "131N", "131C")
  expect_equal(mz_to_tmt_channel(mz, "TMT-11plex"), expected)
})


test_that("mz_to_tmt_channel works for TMTpro-10plex", {
  mz <- c(126.127726, 127.124761, 128.128116, 129.131471, 130.134826,
          131.138181, 132.141536, 133.144891, 134.148246, 135.151601)
  expected <- c("126", "127N", "128N", "129N", "130N",
                "131N", "132N", "133N", "134N", "135N")
  expect_equal(mz_to_tmt_channel(mz, "TMTpro-10plex"), expected)
})


test_that("mz_to_tmt_channel works for TMTpro-16plex", {
  mz <- c(126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
          129.131471, 129.137790, 130.134825, 130.141145, 131.138180,
          131.144499, 132.141534, 132.147853, 133.144888, 133.151207,
          134.148242)
  expected <- c("126", "127N", "127C", "128N", "128C",
                "129N", "129C", "130N", "130C", "131N",
                "131C", "132N", "132C", "133N", "133C",
                "134N")
  expect_equal(mz_to_tmt_channel(mz, "TMTpro-16plex"), expected)
})


test_that("mz_to_tmt_channel works for TMTpro-18plex", {
  mz <- c(126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
          129.131471, 129.137790, 130.134825, 130.141145, 131.138180,
          131.144499, 132.141534, 132.147853, 133.144888, 133.151207,
          134.148242, 134.154562, 135.151597)
  expected <- c("126", "127N", "127C", "128N", "128C",
                "129N", "129C", "130N", "130C", "131N",
                "131C", "132N", "132C", "133N", "133C",
                "134N", "134C", "135N")
  expect_equal(mz_to_tmt_channel(mz, "TMTpro-18plex"), expected)
})
