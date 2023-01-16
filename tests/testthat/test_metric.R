test_that("CONGA works", {
  expect_equal(conga(x=c(100,110,120,90,80), times=c(0,60,120,180,240), n=1, s=1, method="manuscript"), sqrt(1100/3))
})

test_that("LI works", {
  expect_equal(li(x=c(100,110,120,90,80), times=c(0,60,120,180,240), k=60, s=1), 300)
})

test_that("J-index works", {
  expect_equal(j_index(x=c(100,110,120,90,80),"mg"), (100+sqrt(1000/4))^2/1000)
})

test_that("BGI works", {
  expect_equal(bgi(x=c(300,300,300,300,300),"mg")$lbgi, 0)
  expect_equal(bgi(x=c(100,100,100,100,100),"mg")$hbgi, 0)
})

test_that("GRADE works", {
  expect_equal(grade(x=c(100,110,120,90,80))$hypo, 0)
  expect_equal(grade(x=c(100,110,120,90,80))$eu, 1)
  expect_equal(grade(x=c(100,110,120,90,80))$hyper, 0)
})

test_that("MODD works", {
  expect_equal(modd(x=c(100,110,120,90,80), times=seq(0,1440*4,1440)), 15)
})

test_that("MAGE works", {
  expect_equal(mage(x=c(rep(100,10),rep(120,10),105,85), times=seq(0,1260,60)), 20)
})

test_that("MAG works", {
  expect_equal(mag(x=c(100,110,120,90,80), times=c(0,60,120,180,240)), 15)
})

test_that("CV works", {
  expect_equal(cv(x=c(100,100,100,100,100), times=c(0,60,120,180,240)), 0)
})

test_that("SD works", {
  expect_equal(st_dev(x=c(100,110,120,90,80), times=c(0,60,120,180,240)), sd(c(100,110,120,90,80)))
})

test_that("AUC works", {
  expect_equal(cgm_auc(x=c(100,110,120,90,80), times=c(0,60,120,180,240)), 7200)
})

test_that("TIR works", {
  expect_equal(tir(x=c(100,110,120,90,80), low=105, high=125), 40)
})

test_that("num_events works", {
  expect_equal(num_events(x=c(100,110,120,90,80), times=c(0,60,120,180,240), thresh=55, len=60), 0)
})

test_that("GVP works", {
  expect_equal(gvp(x=c(100,100,100,100,100), times=c(0,60,120,180,240)), 0)
})

test_that("DT works", {
  expect_equal(dist_travelled(x=c(100,110,120,90,80)), 60)
})

test_that("time_on works", {
  expect_equal(time_on(times=c(0,5,10,30,35,40), 5), 25)
})
