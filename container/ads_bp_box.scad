// Dimensions in mm
// Use M2.5x12mm screws

con_wid = 5.0;
con_len = 26.0;
con_hgt = 9.0;
brd_hgt = 1.5;

eps = 1.5;                  // Wall width
tol = 0.5;                  // Extruder tolerance
clip_wid = eps;
clip_len = 20.0;
clip_hgt = eps;

inner_wid = 62;
inner_hgt = 45;
inner_len = 130;

module jack_mount() {
    $fn=20;
    rotate([0,-90,0])
        translate([0,0,-2*eps]) cylinder(d=6.5,h=3*eps);
}

module female_clip() {
    clip_tol = eps/2;
    
    translate([-clip_tol/2, -clip_tol/2, 0])
        cube([clip_len+clip_tol,clip_wid+clip_tol,3]);
}

module male_clip() {
    // Rectangular clip
    translate([0,0,-eps]) cube([clip_len, clip_wid, 2*clip_hgt]);
    
    // Taper lower edge into box
    hull() {
        $fn=16;
        translate([0,eps/2,-eps]) 
            rotate([0,90,0]) cylinder(d=eps, h=clip_len);
        translate([0,3*eps/2,-eps-5]) 
            rotate([0,90,0]) cylinder(d=eps, h=clip_len);
    }
}

module screw_mount() {
    // Mounting hole for M2.5x12mm screw
    $fn=20;
    translate([0,0,-12]) {difference() {
        union() {
            hull() {
                // Screw mount outer boundary
                cylinder(d=3*eps, h=12);
                sphere(d=3*eps);
                // Taper bottom edge into box
                translate([-eps-1-1.25,-1.5*eps,-6]) 
                    cube([eps+1,3*eps,18]);
            }
        }
        // Screw hole
        cylinder(d=2.5, h=13);
    } }
}

module outer_box() {
    // Round edges and corners
    module bounding() {
        // This needs to have vectors just slightly larger than eps
        sphere(r=eps+0.013);
    }
    
    $fn=32;
    hull() {
        translate([0,           0,          eps])       bounding();
        translate([0,           inner_wid,  eps])       bounding();
        translate([inner_len,   inner_wid,  eps])       bounding();
        translate([inner_len,   0,          eps])       bounding();
        translate([0,           0,          inner_hgt]) bounding();
        translate([0,           inner_wid,  inner_hgt]) bounding();
        translate([inner_len,   inner_wid,  inner_hgt]) bounding();
        translate([inner_len,   0,          inner_hgt]) bounding();
    }
}

module lower_box() {
    difference() {
        outer_box();
        
        // Inner space
        translate([0,-2*eps,eps]) 
            cube([inner_len,inner_wid+2*eps,inner_hgt+eps]);
        
        // USB holes
        translate([-1.55, 41+eps,   24.5+eps]) cube([1.6,12,10]);  // Comm
        translate([-1.55, 13.5+eps, 18+eps])   cube([1.6,12,10]);  // Dev
        translate([-1.55, 41.5+eps, 6.7+eps])  cube([1.6,12,10]);  // Charge
        
        // Lead holes
        translate([inner_len,13,34]) jack_mount();
        translate([inner_len,25,34]) jack_mount();
        translate([inner_len,37,34]) jack_mount();
        translate([inner_len,49,34]) jack_mount();
        
        translate([inner_len,13,22]) jack_mount();
        translate([inner_len,25,22]) jack_mount();
        translate([inner_len,37,22]) jack_mount();
        translate([inner_len,49,22]) jack_mount();
        
        translate([inner_len,31,10]) jack_mount();  // Reference
        
        // Bottom clip holes
        translate([30,0,-1]) female_clip();
        translate([80,0,-1]) female_clip();
        
        // Shave off edges
        translate([0,inner_wid-eps,inner_hgt-0.5])
            cube([inner_len, 3*eps, 2*eps]);
    }
    
    // Holders for board standoffs
    translate([31.5,4.5,0])
    difference() {
        cube([27+2*eps,6+2*eps,6+eps]);
        translate([eps,eps,eps]) cube([27,6,7]);
    }
    translate([31.5,47.5,0])
    difference() {
        cube([27+2*eps,6+2*eps,6+eps]);
        translate([eps,eps,eps]) cube([27,6,7]);
    }

    // Male clips
    translate([30,inner_wid-eps,inner_hgt]) male_clip();
    translate([80,inner_wid-eps,inner_hgt]) male_clip();
    
    // Screw mounts
    translate([2.5,           1.5*eps+0.5,       inner_hgt-0.5]) 
        screw_mount();
    translate([inner_len-2.5, 1.5*eps+0.5,       inner_hgt-0.5]) 
        rotate([0,0,180]) screw_mount();
    translate([2.5,           inner_wid-1.5*eps, inner_hgt-0.5]) 
        screw_mount();
    translate([inner_len-2.5, inner_wid-1.5*eps, inner_hgt-0.5]) 
        rotate([0,0,180]) screw_mount();
}

module upper_box() {
    difference() {
        outer_box();
        
        // Inner box
        translate([-eps,0,-2*eps]) 
            cube([inner_len+2*eps,inner_wid+2*eps,inner_hgt+2*eps]);
        
        // Shave edges
        translate([inner_len-0.5,-2*eps,-2*eps]) 
            cube([2*eps,inner_wid+eps*4,inner_hgt+eps*4]);
        translate([-2*eps+0.5,-2*eps,-2*eps])
            cube([2*eps,inner_wid+eps*4,inner_hgt+eps*4]);
        translate([-eps,-2*eps,-eps+0.5])         
            cube([inner_len+2*eps,3*eps,2*eps]);
        
        // Power switch hole
        translate([16.5+eps,-2*eps,10+eps]) cube([10,3*eps,8]);
        
        // Female clips
        translate([30,inner_wid-eps,inner_hgt-1]) female_clip();
        translate([80,inner_wid-eps,inner_hgt-1]) female_clip();
        
        // Screw holes
        $fn=20;
        translate([2.5,           1.5*eps+0.5,       inner_hgt-eps]) 
            cylinder(d=2.5,h=3*eps);
        translate([inner_len-2.5, 1.5*eps+0.5,       inner_hgt-eps]) 
            cylinder(d=2.5,h=3*eps);
        translate([2.5,           inner_wid-1.5*eps, inner_hgt-eps]) 
            cylinder(d=2.5,h=3*eps);
        translate([inner_len-2.5, inner_wid-1.5*eps, inner_hgt-eps]) 
            cylinder(d=2.5,h=3*eps);
    }
    
    // Male clips
    translate([30,eps,eps]) rotate([180,0,0]) male_clip();
    translate([80,eps,eps]) rotate([180,0,0]) male_clip();
}


lower_box();
translate([inner_len,inner_wid+8,inner_hgt+eps]) rotate([180,0,180]) 
    upper_box();
