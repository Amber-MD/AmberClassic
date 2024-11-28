#!/usr/bin/perl

# assign voxels to Bragg or neighboring 

<>;<>;  # skip header
$i=0;
while(<>){
   ($h, $k, $l, $qin, $Idiff, $fav, $phiav) = split(' ', $_);
   $Idiff /= 27.;
   $fav2 = ($fav * $fav)/ 27.;

   # find the neighboring bragg indices:
   if( $h % 3 == 0 ) { $hb = $h/3; }
   if( $k % 3 == 0 ) { $kb = $k/3; }
   if( $l % 3 == 0 ) { $lb = $l/3; }
   if( $h % 3 == 1 ) { $hb = ($h - 1)/3; }
   if( $k % 3 == 1 ) { $kb = ($k - 1)/3; }
   if( $l % 3 == 1 ) { $lb = ($l - 1)/3; }
   if( $h % 3 == 2 ) { $hb = ($h + 1)/3; }
   if( $k % 3 == 2 ) { $kb = ($k + 1)/3; }
   if( $l % 3 == 2 ) { $lb = ($l + 1)/3; }

   $key = "$hb\t$kb\t$lb";

   if( $h % 3 == 0 && $k % 3 == 0 && $l % 3 == 0){
       # print "found Bragg: $i   $h  $k  $l  $key\n";
       $diffuse{$key} = $Idiff;
       $bragg{$key} = $fav2;
       $q{$key} = $qin;
   } else {
       $nextto{$key} += $Idiff;
   }
   $i++;
}

print "h\tk\tl\tq\tIdiff\tIbragg\tbetween\tq2B\thalo\n";
print "4N\t4N\t4N\t10N\t12N\t12N\t12N\t12N\t12N\n";

foreach $key (keys %nextto){
    if (exists $q{$key}){
       $nextto{$key} /= 26.;
       print "$key\t$q{$key}\t$diffuse{$key}\t$bragg{$key}\t$nextto{$key}\t";
       $q2b = $q{$key} * $q{$key} * $bragg{$key};
       $diff = $diffuse{$key} - $nextto{$key};
       print "$q2b\t$diff\n";
    }
}

