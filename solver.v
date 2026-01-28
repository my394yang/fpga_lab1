//========================================================================
// Hardware Layout
//========================================================================


module Data_Path
(
  (* keep=1 *) input  logic [26:0] x_0,
  (* keep=1 *) input  logic [26:0] y_0,
  (* keep=1 *) input  logic [26:0] z_0,
  (* keep=1 *) input  logic [26:0] rho,
  (* keep=1 *) input  logic [26:0] beta,
  (* keep=1 *) input  logic [26:0] sigma
  //(* keep=1 *) input  logic [26:0] dt i think will be given as input


  (* keep=1 *) output  logic [26:0] xk,
  (* keep=1 *) output  logic [26:0] yk,
  (* keep=1 *) output  logic [26:0] zk,

);
//confused on the integrator and how this plays into it, sorry 

/*
module integrator(out,funct,InitialOut,clk,reset);
	output signed [26:0] out; 		//the state variable V
	input signed [26:0] funct;      //the dV/dt function
	input clk, reset;
	input signed [26:0] InitialOut;  //the initial state variable V
	
	wire signed	[26:0] out, v1new ;
	reg signed	[26:0] v1 ;
	
	always @ (posedge clk) 
	begin
		if (reset==0) //reset	
			v1 <= InitialOut ; // 
		else 
			v1 <= v1new ;	
	end
	assign v1new = v1 + funct ;
	assign out = v1 ;
endmodule
*/

module Adder_27b_RTL
(
  (* keep=1 *) input  logic [26:0] in0,
  (* keep=1 *) input  logic [26:0] in1,
  (* keep=1 *) output logic [26:0] sum
);

assign sum = in0 + in1;
  
endmodule

module Sub_27b_RTL
(
  (* keep=1 *) input  logic [26:0] in0,
  (* keep=1 *) input  logic [26:0] in1,
  (* keep=1 *) output logic [26:0] diff
);

assign diff = in0 - in1;
  
endmodule

module Mult_27x27b_RTL
(
  (* keep=1 *) input  logic [26:0] in0,
  (* keep=1 *) input  logic [26:0] in1,
  (* keep=1 *) output logic [26:0] prod
);

  assign prod = ( in0*in1 );

endmodule
/*
module signed_mult (out, a, b);
	output 	signed  [26:0]	out;
	input 	signed	[26:0] 	a;
	input 	signed	[26:0] 	b;
	// intermediate full bit length
	wire 	signed	[53:0]	mult_out;
	assign mult_out = a * b;
	// select bits for 7.20 fixed point
	assign out = {mult_out[53], mult_out[45:20]};
endmodule
*/

module Reg_RTL
(
  (* keep=1 *) input  logic clk,
  (* keep=1 *) input  logic [26:0] d,
  (* keep=1 *) input  logic en,
  (* keep=1 *) input  logic rst,
  (* keep=1 *) output logic [26:0] q
);

  always_ff @( posedge clk ) begin

      if ( rst )
        q <= 27'b0;
      else if ( en )
        q <= d;

end

endmodule

//mux
module mux_RTL
(
  (* keep=1 *) input  logic [26:0] in0,
  (* keep=1 *) input  logic [26:0] in1,
  (* keep=1 *) input  logic  rst
  (* keep=1 *) output logic  conditions

);

assign conditions = rst ? in1 : in0;

endmodule


//dx/dt= sigma*y-x
logic [26:0] xy_diff;
logic [26:0] dx_dt;
logic [26:0] dx;
logic [26:0] dt;
logic [26:0] d_reg_xk;
logic [26:0] mux_in0_x;
logic rst; //needs to be controlled
//logic [26:0] x_0; //needs to be controlled
logic [26:0] d_reg_x;
logic [4:0] clk; // do we need a clock divider to set clk???
/*

// clock divider to slow system down for testing
reg [4:0] count;
// analog update divided clock
always @ (posedge CLOCK_50) 
begin
        count <= count + 1; 
end	
assign AnalogClock = (count==0);
*/
//assign dt= 1/frequency; something like this 

Sub_27b_RTL dx_calc 
(
  .in0(yk),
  .in1(xk),
  .diff(xy_diff)
);

Mult_27bx27b_RTL dx_dt
(
  .in0(sigma),
  .in1(xy_diff),
  .prod(dx_dt)
);

Mult_27bx27b_RTL dx_calc
(
  .in0(dx_dt),
  .in1(dt),
  .prod(dx)
);

Adder_27b_RTL mux_input_x
(
  .in0(xk),
  .in1(dx),
  .sum(mux_in0_x)
);

Mux_27b_RTL reg_input_x
(
  .in0(mux_in0),
  .in1(x_0),
  .rst(rst),
  .conditions(d_reg_x)
);

Reg_27b_RTL xk_next
(
  .clk(clk)
  .d(d_reg_x),
  .q(xk) //needs to be put into xk
);
//assign xy_diff= yk-xk; // adder
//assign dx_dt= sigma* xy_diff; //multiplier
//assign dx= dx_dt*dt; //multiplier


//yk
//dx/dt= sigma*y-x
logic [26:0] pz_diff;
logic [26:0] dy_dt;
logic [26:0] dy_dt_y;
logic [26:0] dy;
logic [26:0] d_reg_yk;
logic [26:0] mux_in0_y;
logic rst; //needs to be controlled
logic [26:0] d_reg_y;
/*

// clock divider to slow system down for testing
reg [4:0] count;
// analog update divided clock
always @ (posedge CLOCK_50) 
begin
        count <= count + 1; 
end	
assign AnalogClock = (count==0);
*/
//assign dt= 1/frequency; something like this 

Sub_27b_RTL dy_calc 
(
  .in0(rho),
  .in1(zk),
  .diff(pz_diff)
);

Mult_27bx27b_RTL dy_dt_plus_y
(
  .in0(xk),
  .in1(pz_diff),
  .prod(dy_dt_y) //is dy/dt+y
);

Sub_27b_RTL dy_dt 
(
  .in0(dy_dt_y),
  .in1(yk),
  .diff(dy_dt) //is dy/dt
);

Mult_27bx27b_RTL dy_calc
(
  .in0(dy_dt),
  .in1(dt),
  .prod(dy)
);

Adder_27b_RTL mux_input_y
(
  .in0(yk),
  .in1(dy),
  .sum(mux_in0_y)
);

Mux_27b_RTL reg_input_y
(
  .in0(mux_in0_y),
  .in1(x_0),
  .rst(rst),
  .conditions(d_reg_x)
);

Reg_27b_RTL yk_next
(
  .clk(clk)
  .d(d_reg_y),
  .q(yk) //needs to be put into yk
);


//zk
logic [26:0] xy_prod;
logic [26:0] dz_dt;
logic [26:0] beta_z_prod;
logic [26:0] dz;
logic [26:0] d_reg_zk;
logic [26:0] mux_in0_z;
logic rst; //needs to be controlled
logic [26:0] d_reg_z;

Mult_27bx27b_RTL xy
(
  .in0(xk),
  .in1(yk),
  .prod(xy_prod) 
);

Mult_27bx27b_RTL Beta_z
(
  .in0(beta),
  .in1(zk),
  .prod(beta_z_prod)
);

Sub_27b_RTL dz_dt 
(
  .in0(xy_prod),
  .in1(beta_z_prod),
  .diff(dz_dt) //is dz/dt
);

Mult_27bx27b_RTL dz_calc
(
  .in0(dt),
  .in1(dz_dt),
  .prod(dz)
);

Adder_27b_RTL mux_input_z
(
  .in0(zk),
  .in1(dz),
  .sum(mux_in0_z)
);

Mux_27b_RTL reg_input_z
(
  .in0(mux_in0_z),
  .in1(z_0),
  .rst(rst),
  .conditions(d_reg_z)
);

Reg_27b_RTL zk_next
(
  .clk(clk)
  .d(d_reg_z),
  .q(zk) //needs to be put into zk
);


endmodule

`endif

