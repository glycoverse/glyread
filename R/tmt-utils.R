mz_to_tmt_channel <- function(x, tmt_type) {
  tmt_type <- rlang::arg_match(
    tmt_type,
    c("TMT-6plex", "TMT-10plex", "TMT-11plex",
      "TMTpro-10plex", "TMTpro-16plex", "TMTpro-18plex")
  )

  tmt11_channels <- c(
    "126" = 126.1277,
    "127N" = 127.1248,
    "127C" = 127.1311,
    "128N" = 128.1281,
    "128C" = 128.1344,
    "129N" = 129.1315,
    "129C" = 129.1378,
    "130N" = 130.1348,
    "130C" = 130.1411,
    "131N" = 131.1382,
    "131C" = 131.1445
  )
  tmtpro18_channels <- c(
    "126" = 126.1277,
    "127N" = 127.1248,
    "127C" = 127.1311,
    "128N" = 128.1281,
    "128C" = 128.1344,
    "129N" = 129.1315,
    "129C" = 129.1378,
    "130N" = 130.1348,
    "130C" = 130.1411,
    "131N" = 131.1382,
    "131C" = 131.1445,
    "132N" = 132.1415,
    "132C" = 132.1479,
    "133N" = 133.1449,
    "133C" = 133.1512,
    "134N" = 134.1482,
    "134C" = 134.1546,
    "135N" = 135.1516
  )

  if (tmt_type == "TMT-6plex") {
    channels <- tmt11_channels[c("126", "127N", "128C", "129N", "130C", "131N")]
    names(channels) <- c("126", "127", "128", "129", "130", "131")
  } else if (tmt_type == "TMT-10plex") {
    channels <- tmt11_channels[1:10]
    names(channels)[[10]] <- "131"
  } else if (tmt_type == "TMT-11plex") {
    channels <- tmt11_channels
  } else if (tmt_type == "TMTpro-10plex") {
    channels <- tmtpro18_channels[c("126", paste0(127:135, "N"))]
  } else if (tmt_type == "TMTpro-16plex") {
    channels <- tmtpro18_channels[1:16]
  } else if (tmt_type == "TMTpro-18plex") {
    channels <- tmtpro18_channels
  }

  x <- round(as.double(x), 4)
  names(channels)[match(x, channels)]
}
