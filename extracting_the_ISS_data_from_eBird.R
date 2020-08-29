### extracting the ISS data from eBird

library(auk)

# untar("E:/eBird_all/ebd_sampling_relJun-2020.tar",exdir = "E:/eBird_all")
# untar("E:/eBird_all/ebd_relJun-2020.tar",exdir = "E:/eBird_all")
### setting path to eBird full dataset
auk::auk_set_ebd_path("E:/eBird_all",overwrite = T)


ebd <- auk_ebd("E:/eBird_all/ebd_relJun-2020.txt", 
               file_sampling = "E:/eBird_all/ebd_sampling_relJun-2020.txt")



ebd_filters <- ebd %>% 
  # auk_species("Wood Thrush") %>% 
  # # southeastern coastal plain bcr
  # auk_bcr(bcr = 27) %>% 
  # fall surveys only
  auk_date(date = c("*-07-01", "*-12-31")) %>% 
  # restrict to the ISS protocol
  auk_protocol(protocol = c("International Shorebird Survey (ISS)")) 



data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
f_ebd <- file.path(data_dir, "ebd_ISS.txt")
f_sampling <- file.path(data_dir, "ebd_ISS.txt")

# only run if the files don't already exist
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
}


