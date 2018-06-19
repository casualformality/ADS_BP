// Dimensions in mm
// Use M2.5x12mm screws

con_wid = 5.0;
con_len = 26.0;
con_hgt = 9.0;
brd_hgt = 1.5;

eps = 1.5;
tol = eps/2;
clip_wid = eps;
clip_len = 20.0;
clip_hgt = eps;

module battery_board() {
    brd_wid = 59.0;
    brd_len = 66.0;
    
    color("green") {
        translate([29,5,0])
        cube([con_len, con_wid, con_hgt]);

        translate([29,48,0])
        cube([con_len, con_wid, con_hgt]);

        translate([0,0,9])
        cube([brd_len, brd_wid, brd_hgt]);
    }
    
    color("silver") {
        // usb port
        translate([-1,43,10.5])
        cube([6,9,3]);
        
        // power switch
        translate([14.25,-11,10.5])
        cube([9,21,6]);
    }
}

module tiva_board() {
    brd_wid = 51.0;
    brd_len = 66.0;
    
    color("red") {
        translate([29,1,0])
        cube([con_len, con_wid, con_hgt]);

        translate([29,44,0])
        cube([con_len, con_wid, con_hgt]);

        translate([0,0,9])
        cube([brd_len, brd_wid, brd_hgt]);
    }
    
    color("silver") {
        // usb port
        translate([-1,11,10.5])
        cube([6,9,3]);
    }
}

module comm_board() {
    brd_wid = 52.0;
    brd_len = 66.0;
        
    color("violet") {
        translate([22,1.5,0])
        cube([con_len, con_wid, con_hgt]);

        translate([22,44.5,0])
        cube([con_len, con_wid, con_hgt]);
        
        translate([22,1.5,10.5])
        cube([con_len, con_wid, con_hgt]);

        translate([22,44.5,10.5])
        cube([con_len, con_wid, con_hgt]);
        
        translate([0,0,9])
        cube([brd_len, brd_wid, brd_hgt]); 
    }
    
    color("silver") {
        // input ports
        translate([59,14,10.5])
        cube([14,21,5]);
        
        // usb port
        // hole should be 12x10mm
        translate([-1,38.5,5])
        cube([10,8,4]);
    }
}

module board_stack() {
    battery_board();
    translate([0,4,12]) tiva_board();
    translate([7,3.5,23]) comm_board();
}

module jack_mount() {
    $fn=20;
    rotate([0,-90,0])
    union() {
        // 1mm diameter tolerance for circular pieces
        //translate([0,0,-4.1]) cylinder(r=4.5,h=2.1);
        translate([0,0,-4]) cylinder(r=3.5,h=4.1);
    }
}

module female_clip() {
    translate([-tol/2, -tol/2, 0])
        cube([clip_len+tol,clip_wid+tol,3]);
}

module male_clip() {
    translate([0,0,-eps]) cube([clip_len, clip_wid, 2*clip_hgt]);
    hull() {
        $fn=16;
        translate([0,eps/2,-eps]) 
            rotate([0,90,0]) cylinder(d=eps, h=clip_len);
        translate([0,3*eps/2,-eps-5]) 
            rotate([0,90,0]) cylinder(d=eps, h=clip_len);
    }
}

// Mounting hole for M2.5x12mm screw
module screw_mount() {
    $fn=20;
    difference() {
        union() {
            hull() {
                cylinder(d=3*eps, h=12);
                sphere(d=3*eps);
                translate([-3,-1.5*eps,-5-1.5*eps]) 
                    cube([1.5*eps,3*eps,19.5]);
            }
        }
        cylinder(d=2.5, h=13);
    }
}

inner_wid = 62;
inner_hgt = 47;
inner_len = 120;
outer_len = inner_len+eps;

module outer_box() {
    module bounding() {
        rotate([90,0,0]) sphere(r=eps*1.01);
    }
    
    $fn=32;
    hull() {
        translate([0,0,eps]) bounding();
        translate([0,inner_wid,eps]) bounding();
        translate([outer_len-eps,inner_wid,eps]) bounding();
        translate([outer_len-eps,0,eps]) bounding();
        translate([0,0,inner_hgt]) bounding();
        translate([0,inner_wid,inner_hgt]) bounding();
        translate([outer_len-eps,inner_wid,inner_hgt]) bounding();
        translate([outer_len-eps,0,inner_hgt]) bounding();
    }
}

module lower_box() {
    difference() {
        outer_box();
        
        // Inner space
        translate([0,-2,eps]) cube([inner_len,64,47+eps]);
        
        // USB holes
        translate([-1.55,40+eps,25.5+eps]) cube([1.6,12,10]);
        translate([-1.55,13.5+eps,19+eps]) cube([1.6,12,10]);
        translate([-1.55,41.5+eps,7.7+eps]) cube([1.6,12,10]);
        
        // Power switch hole
        translate([14.5+eps,-1.55,10+eps]) cube([10,1.6,8]);
        
        // Lead holes
        translate([inner_len,13,34]) jack_mount();
        translate([inner_len,25,34]) jack_mount();
        translate([inner_len,37,34]) jack_mount();
        translate([inner_len,49,34]) jack_mount();
        translate([inner_len,13,22]) jack_mount();
        translate([inner_len,25,22]) jack_mount();
        translate([inner_len,37,22]) jack_mount();
        translate([inner_len,49,22]) jack_mount();
        translate([inner_len,13,10]) jack_mount();
        
        // Bottom clip holes
        translate([20,0,-1]) female_clip();
        translate([70,0,-1]) female_clip();
        
        // Shave off edges
        translate([0,inner_wid-eps/2,46.75]) 
            cube([inner_len,2*eps,2*eps]);
    }
    
    // Holders for standoffs
    translate([31.5,4.5,0])
    difference() {
        cube([26.5+2*eps,6+2*eps,5+eps]);
        translate([eps,eps,eps]) cube([26.5,6,6]);
    }
    translate([31.5,47.5,0])
    difference() {
        cube([27+2*eps,6+2*eps,5+eps]);
        translate([eps,eps,eps]) cube([27,6,6]);
    }

    // Male clips
    translate([20,inner_wid-eps,inner_hgt]) male_clip();
    translate([70,inner_wid-eps,inner_hgt]) male_clip();
    
    // Screw mounts
    translate([2.5,1.5*eps,inner_hgt-12.25]) screw_mount();
    translate([inner_len-2.5,1.5*eps,inner_hgt-12.25]) 
        rotate([0,0,180]) screw_mount();
}

module upper_box() {
    difference() {
        outer_box();
        
        // Inner box
        translate([-2*eps,0,-2*eps]) cube([outer_len+2*eps,64,inner_hgt+2*eps]);
        
        // Shave edges
        translate([inner_len-tol/2,-2*eps,0]) cube([4+eps+tol/2,68,51]);
        translate([-eps-tol/2,-2*eps,0]) cube([eps+tol,68,51]);
        translate([-2*eps,-2*eps,-eps]) cube([outer_len+4*eps,3*eps,tol/2+2*eps]);
        
        // Power switch hole
        translate([16.5+eps,-1.55,9.5+2]) cube([10,1.6,8]);
        
        // Female clips
        translate([20,inner_wid-eps,inner_hgt-1]) female_clip();
        translate([70,inner_wid-eps,inner_hgt-1]) female_clip();
        
        // Screw holes
        $fn=20;
        translate([2.5,1.5*eps,inner_hgt-eps]) cylinder(d=2.5,h=3*eps);
        translate([inner_len-2.5,1.5*eps,inner_hgt-eps]) 
            cylinder(d=2.5,h=3*eps);
    }
    
    // Male clips
    translate([20,eps,eps]) rotate([180,0,0]) male_clip();
    translate([70,eps,eps]) rotate([180,0,0]) male_clip();
}


//translate([eps+2.75,eps,2]) board_stack();
lower_box();
translate([inner_len,inner_wid+8,inner_hgt+eps]) rotate([180,0,180]) 
    upper_box();

