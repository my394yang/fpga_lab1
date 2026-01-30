`timescale 1ns/1ns

module testbench();
	
	reg clk_50, reset;
	
	reg [31:0] index;
	wire signed [26:0] x_out, y_out, z_out;
	
	//Initialize clocks and index
	initial begin
		clk_50 = 1'b0;
		index  = 32'd0;
	end
	
	//Toggle the clock
	always begin
		#10
		clk_50  = !clk_50;
	end
	
	//Initialize and drive signals
	initial begin
		reset  = 1'b0;
		#40 
		reset  = 1'b1;
	end
	
	//Increment index
	always @ (posedge clk_50) begin
		index  <= index + 32'd1;
	end

	//Instantiation of Device Under Test
	lorenz_system DUT (
		.clock(clk_50), 
		.reset(reset),
		.x(x_out),
		.y(y_out),
		.z(z_out)
	);
	
endmodule

/////////////////////////////////////////////////
////////////	Direct Digital Synth	//////////////
//////////////////////////////////////////////////

module lorenz_system(clock, reset, x, y, z);
input clock, reset;
output wire signed [26:0] x, y, z;

// Declare wires for initial conditions and constants:
wire signed [26:0] x_init, y_init, z_init;
wire signed [26:0] sigma, rho, beta;
wire signed [26:0] dt;

// Initial condition integers of -1, 0.1 and 25 converted to 7.20 fixed point convention
assign x_init = -27'sb0000001_00000000000000000000;
assign y_init =  27'sb0000000_00011001100110011001;
assign z_init =  27'sb0011001_00000000000000000000;
assign dt = 27'sb0000000_00000001000000000000;

// Constants
assign sigma = 27'sb0001010_00000000000000000000;
assign rho   = 27'sb0011100_00000000000000000000;
assign beta  = 27'sb0000010_10101010101010101011;

// Increments (outputs from derivative modules)
wire signed [26:0] x_inc, y_inc, z_inc;

// clock divider to slow system down for testing
reg [4:0] count;
wire AnalogClock;

// analog update divided clock
always @ (posedge clock) 
begin
	if (reset == 0) //reset active
		count <= 5'd0;
	else
		count <= count + 1; 
end 

assign AnalogClock = (count==0);

/////////////////////////////////////////////////
// Instantiate derivative computation modules
/////////////////////////////////////////////////

increment_x compute_dx (
	.x_increment_result(x_inc),
	.x(x),
	.y(y),
	.sigma(sigma),
	.dt(dt)
);

increment_y compute_dy (
	.y_increment_result(y_inc),
	.x(x),
	.y(y),
	.z(z),
	.rho(rho),
	.dt(dt)
);

increment_z compute_dz (
	.z_increment_result(z_inc),
	.x(x),
	.y(y),
	.z(z),
	.beta(beta),
	.dt(dt)
);


/////////////////////////////////////////////////
// Instantiate integrators
/////////////////////////////////////////////////

integrator integrate_x (
	.out(x),
	.funct(x_inc),
	.InitialOut(x_init),
	.clk(AnalogClock),
	.reset(reset)
);

integrator integrate_y (
	.out(y),
	.funct(y_inc),
	.InitialOut(y_init),
	.clk(AnalogClock),
	.reset(reset)
);

integrator integrate_z (
	.out(z),
	.funct(z_inc),
	.InitialOut(z_init),
	.clk(AnalogClock),
	.reset(reset)
);

endmodule


/////////////////////////////////////////////////
//// integrator /////////////////////////////////
/////////////////////////////////////////////////

module integrator(out, funct, InitialOut, clk, reset);
output signed [26:0] out;       //the state variable V
input signed [26:0] funct;      //the dV/dt function
input clk, reset;
input signed [26:0] InitialOut;  //the initial state variable V
wire signed [26:0] out, v1new;
reg signed [26:0] v1;

always @ (posedge clk) 
begin
	if (reset==0) //reset   
		v1 <= InitialOut; // 
	else 
		v1 <= v1new;   
end

assign v1new = v1 + funct;
assign out = v1;
endmodule

/////////////////////////////////////////////////
//// signed multiplier of 7.20 format 2'comp/////
/////////////////////////////////////////////////

module signed_mult (out, a, b);
output signed [26:0] out;
input signed [26:0] a;
input signed [26:0] b;
// intermediate full bit length
wire signed [53:0] mult_out;
assign mult_out = a * b;
// select bits for 7.20 fixed point
assign out = {mult_out[53], mult_out[45:20]};
endmodule

/////////////////////////////////////////////////
//// Compute dx/dt///////////////////////////////
/////////////////////////////////////////////////

module increment_x(x_increment_result, x, y, sigma, dt);
// Stores the value of dt dot dx/dt
output signed [26:0] x_increment_result;
input signed [26:0] x, y, sigma, dt;
wire signed [26:0] xy_difference, dt_scaled_diff, dx_dt;

assign xy_difference = y - x;

// Multiplies (x-y) with dt first to avoid overflow error
signed_mult mult_dt_first (
    .out(dt_scaled_diff),
    .a(dt),
    .b(xy_difference)
);

signed_mult mult_sigma (
    .out(x_increment_result),
    .a(sigma),
    .b(dt_scaled_diff)
);
endmodule

/////////////////////////////////////////////////
//// Compute dy/dt///////////////////////////////
/////////////////////////////////////////////////

module increment_y(y_increment_result, x, y, z, rho, dt);
// Stores the value of dt dot dy/dt
output signed [26:0] y_increment_result;
input signed [26:0] x, y, z, rho, dt;
wire signed [26:0] rho_minus_z_scaled, x_times_scaled_term;

// Scales (rho - z) by dt first
assign rho_minus_z_scaled = (rho - z) >>> 8;

// Computes x * dt * (rho - z)
signed_mult mult_x_term (
    .out(x_times_scaled_term),
    .a(x),
    .b(rho_minus_z_scaled)
);

assign y_increment_result = x_times_scaled_term - (y >>> 8);

endmodule

/////////////////////////////////////////////////
//// Compute dz/dt///////////////////////////////
/////////////////////////////////////////////////

module increment_z(z_increment_result, x, y, z, beta, dt);
// Stores the value of dt dot dz/dt
output signed [26:0] z_increment_result;
input signed [26:0] x, y, z, beta, dt;
wire signed [26:0] x_scaled, xy_scaled_result, beta_z_scaled;

// Scales x by dt first to prevent overflow
assign x_scaled = x >>> 8;

// Computes dt*x*y
signed_mult mult_xy_scaled (
    .out(xy_scaled_result),
    .a(x_scaled),
    .b(y)
);

// Computes dt*beta*z
signed_mult mult_beta_z (
    .out(beta_z_scaled),
    .a(beta >>> 8),
    .b(z)
);

assign z_increment_result = xy_scaled_result - beta_z_scaled;

endmodule

