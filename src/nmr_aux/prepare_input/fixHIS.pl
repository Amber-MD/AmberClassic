#!/usr/bin/env perl -n
#
#   take output from reduce, change HIS -> HID, HIE or HIP based on FlipMemo
#
#   first, read in a pdb file line and unpack
#

if( $_ =~ /no HE2:/){ $hidno[ $nhid++ ] = substr( $_, 21, 3);}
if( $_ =~ /no HD1:/){ $hieno[ $nhie++ ] = substr( $_, 21, 3);}
if( $_ =~ /\+bothHN/){ $hipno[ $nhip++ ] = substr( $_, 21, 3);}

if( $_ =~ /^ATOM|^HETATM/){

( $atno, $atname, $alt, $resname, $chainId, $resno, $iCode,
      $x, $y, $z, $occ, $bfact, $element, $charge ) =
    unpack("x7 a4 x a4 a a3 x a a4 a x3 a8 a8 a8 a6 a6 x10 a2 a2",$_);

#
#  do modifications necessary here:
#
$charge="  ";   #reduce doesnt output charge
for( $i=0; $i<$nhid; $i++ ){ $resname = "HID" if $hidno[$i] == $resno;}
for( $i=0; $i<$nhie; $i++ ){ $resname = "HIE" if $hieno[$i] == $resno;}
for( $i=0; $i<$nhip; $i++ ){ $resname = "HIP" if $hipno[$i] == $resno;}
$resno -= 3;
#
#  write back out
#
printf 
"ATOM   %4s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
      $atno,$atname,$alt, $resname,$chainId, $resno, $iCode,
        $x,$y,$z, $occ, $bfact, $element, $charge;

} else {
   print;
}
