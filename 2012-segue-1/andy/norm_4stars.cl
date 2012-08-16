noao
onedspec
obsutil
imred
ccdr
ech
!rm -f red.fits
!rm -f red_multic.fits
!rm -f red_multicim.fits
!rm -f red_multicom.fits
!rm -f red_multicimcom.fits
!rm -f endred.fits
!rm -f a.fits
!rm -f b.fits
!rm -f c.fits
!rm -f d.fits
!rm -f e.fits
!rm -f f.fits
scopy (input="seg11red_multi_ut100318_ut100319.fits",  output="a.fits", aperture="", bands=2, beams="1", format="multispec")
scopy (input="seg11red_multi_ut100321.fits",  output="c.fits", aperture="", bands=2, beams="1", format="multispec")
scopy (input="segue1-11red_multi_ut100307_ut100308.fits",  output="d.fits", aperture="", bands=2, beams="1", format="multispec")
onedspec
scombine ("a.fits,c.fits,d.fits","red.fits",aperture="*",combine="sum",group="aper", reject="ccdclip",lthresh=0.,hthresh=15000.,nlow=0.,rdnoise="ENOISE",gain="EGAIN")
cont (input="red.fits", output="red_multic.fits", order=6, overrid=yes, markrej=no)

imarith (operand1="red.fits", op="/", operand2="red_multic.fits", result="red_multicim.fits")
noao
onedspec
scombine ("red.fits", "red_multicom.fits", combine="sum", group="images")
scombine ("red_multicim.fits", "red_multicimcom.fits",  combine="sum", group="images")
imarith  (operand1="red_multicom.fits", op="/", operand2="red_multicimcom.fits", result="endred.fits")
copy ("endred.fits", "seg11red.fits")
wspectext     (input="seg11red.fits", output="seg11red.txt", header=no)

------------------------------------------

noao
onedspec
obsutil
imred
ccdr
ech
!rm -f blue.fits
!rm -f blue_multic.fits
!rm -f blue_multicim.fits
!rm -f blue_multicom.fits
!rm -f blue_multicimcom.fits
!rm -f endblue.fits
!rm -f a.fits
!rm -f b.fits
!rm -f c.fits
!rm -f d.fits
!rm -f e.fits
!rm -f f.fits
scopy (input="seg11blue_multi_ut100318_ut100319.fits",  output="a.fits", aperture="73-106", bands=2, beams="1", format="multispec")
scopy (input="seg11blue_multi_ut100321.fits",  output="c.fits", aperture="73-106", bands=2, beams="1", format="multispec")
scopy (input="segue1-11blue_multi_ut100307_ut100308.fits",  output="d.fits", aperture="73-106", bands=2, beams="1", format="multispec")
onedspec
scombine ("a.fits,c.fits,d.fits","blue.fits",aperture="*",combine="sum",group="aper", reject="ccdclip",lthresh=0.,hthresh=15000.,nlow=0.,rdnoise="ENOISE",gain="EGAIN")
cont (input="blue.fits", output="blue_multic.fits", order=6, overrid=yes, markrej=no)

imarith (operand1="blue.fits", op="/", operand2="blue_multic.fits", result="blue_multicim.fits")
noao
onedspec
scombine ("blue.fits", "blue_multicom.fits", combine="sum", group="images")
scombine ("blue_multicim.fits", "blue_multicimcom.fits",  combine="sum", group="images")
imarith  (operand1="blue_multicom.fits", op="/", operand2="blue_multicimcom.fits", result="endblue.fits")
copy ("endblue.fits", "seg11blue.fits")
wspectext      (input="seg11blue.fits", output="seg11blue.txt", header=no)

=====================================================

NOT USED! too few counts!
scopy (input="seg11blue_multi_ut100322_ut100323.fits",  output="b.fits", aperture="", bands=2, beams="1", format="multispec")
scopy (input="seg11red_multi_ut100322_ut100323.fits",  output="b.fits", aperture="", bands=2, beams="1", format="multispec")


------------------------
tests
splot ("seg11blue_multi_ut100318_ut100319.fits",line=1, band=2)
splot ("seg11blue_multi_ut100321.fits",line=2, band=2)
splot ("segue1-11blue_multi_ut100307_ut100308.fits",line=3, band=2)

!rm -f foo.fits
scopy (input="segue1-11blue_multi_ut100307_ut100308.fits",  output="foo.fits", aperture="73-106", bands="2", beams="1", format="multispec")
!rm -f foo2.fits
scopy (input="foo.fits",  output="foo2.fits", aperture="73-106", bands="", beams="", format="multispec")
splot (foo.fits,line=3, band=2)
