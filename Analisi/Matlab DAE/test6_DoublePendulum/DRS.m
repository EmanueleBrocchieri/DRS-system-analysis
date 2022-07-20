%
% Matlab code for the Course:
%
%     Modelling and Simulation of Mechatronics System
%
% by
% Enrico Bertolazzi
% Dipartimento di Ingegneria Industriale
% Universita` degli Studi di Trento
% email: enrico.bertolazzi@unitn.it
%
classdef DRS < DAC_ODEclass
  properties (SetAccess = protected, Hidden = true)
    L0;
    L1;
    L2;
    L4;
    L5;
    xP1;
    yP1;
    xP5;
    yP5;
    alpha;
    L_cylinder;
    m_cylinder;
    Iz_cylinder;
    L_piston;
    m_piston;
    Iz_piston;
    x_triangle;
    y_triangle;
    m_triangle;
    Iz_triangle;
    x_flap;
    y_flap;
    m_flap;
    Iz_flap;
    gravity;
    % baumgarte stabilization
    omega;
    eta;
    npos;
  end
  methods (Access = private)
    function res__bigRHS = bigRHS( self, t, vars__ )
      % extract parameters
      L0 = self.L0;
      L1  = self.L1;
      L2  = self.L2;
      L4  = self.L4;
      L5  = self.L5;
      xP1  = self.xP1;
      yP1  = self.yP1;
      xP5  = self.xP5;
      yP5  = self.yP5;
      alpha  = self.alpha;
      L_cylinder  = self.L_cylinder;
      m_cylinder  = self.m_cylinder;
      Iz_cylinder  = self.Iz_cylinder;
      L_piston  = self.L_piston;
      m_piston  = self.m_piston;
      Iz_piston  = self.Iz_piston;
      x_triangle  = self.x_triangle;
      y_triangle  = self.y_triangle;
      m_triangle  = self.m_triangle;
      Iz_triangle  = self.Iz_triangle;
      x_flap  = self.x_flap;
      y_flap  = self.y_flap;
      m_flap  = self.m_flap;
      Iz_flap  = self.Iz_flap;
      g  = self.gravity;
      omega  = self.omega;
      eta  = self.eta;
      npos  = self.npos;
    

      % extract states
      x1 = vars__(1);
      y1 = vars__(2);
      x2 = vars__(3);
      y2 = vars__(4);
      u1 = vars__(5);
      v1 = vars__(6);
      u2 = vars__(7);
      v2 = vars__(8);

      % evaluate function
      res__2 = -m1 * g;
      res__4 = -m2 * g;
      t3 = L1 ^ 2;
      t4 = x1 ^ 2;
      t5 = y1 ^ 2;
      t7 = omega ^ 2;
      t15 = u1 ^ 2;
      t16 = 2 * t15;
      t17 = v1 ^ 2;
      res__5 = t7 * (t3 - t4 - t5) - 4 * eta * (x1 * u1 + y1 * v1) * omega - t16 - 2 * t17;
      t21 = x2 ^ 2;
      t27 = -x2 + x1;
      t32 = v1 - v2;
      t40 = u2 ^ 2;
      t42 = t32 ^ 2;
      res__6 = t7 * (-t4 + 2 * x2 * x1 - t21 + (L2 + y1 - y2) * (L2 - y1 + y2)) - 4 * omega * eta * (u1 * t27 - u2 * t27 + t32 * (-y2 + y1)) - t16 + 4 * u2 * u1 - 2 * t40 - 2 * t42;

      % store on output
      res__bigRHS = zeros(6,1);
      res__bigRHS(2) = res__2;
      res__bigRHS(4) = res__4;
      res__bigRHS(5) = res__5;
      res__bigRHS(6) = res__6;

    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res__JbigRHS = JbigRHS( self, t, vars__ )
      % extract parameters
      g     = self.gravity;
      m1    = self.m1;
      m2    = self.m2;
      L1    = self.L1;
      L2    = self.L2;
      omega = self.omega;
      eta   = self.eta;

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta2 = vars__(3);
      theta3 = vars__(4);
      theta4 = vars__(5);
      theta5 = vars__(6);
      s__dot = vars__(7);
      theta1__dot = vars__(8);
      theta3__dot = vars__(9);
      theta4__dot = vars__(10);
      theta5__dot = vars__(11);
      theta2__dot = vars__(12);

      % evaluate function
      t1 = theta1__dot ^ 2;
      res__1_1 = 5 * t1;
      t3 = 10 * s + 1602;
      res__1_8 = theta1__dot * t3;
      res__2_1 = -10 * s__dot * theta1__dot;
      t6 = -t3;
      res__2_7 = theta1__dot * t6;
      res__2_8 = s__dot * t6;
      t7 = sin(theta5);
      t8 = theta5__dot ^ 2;
      t9 = t8 * t7;
      t11 = theta4__dot ^ 2;
      t12 = cos(theta4);
      t13 = t12 * t11;
      t15 = cos(theta5);
      t16 = t15 * t8;
      t19 = cos(theta3);
      t21 = sin(theta3);
      t22 = sin(theta4);
      t23 = t11 * t22;
      res__4_4 = t19 * (-840 * t9 - 89600 * t13 - 125580 * t16) - 89600 * (t23 + 0.897e3 / 0.640e3 * t9 - 0.3e1 / 0.320e3 * t16) * t21;
      t31 = t12 * t19 + t22 * t21;
      res__4_5 = 89600 * t11 * t31;
      t34 = t15 + 0.2e1 / 0.299e3 * t7;
      t37 = t15 - 0.299e3 / 0.2e1 * t7;
      t40 = t19 * t34 - 0.2e1 / 0.299e3 * t37 * t21;
      res__4_6 = 125580 * t8 * t40;
      t44 = -t21 * t12 + t22 * t19;
      res__4_10 = 179200 * theta4__dot * t44;
      t49 = t19 * t37 + 0.299e3 / 0.2e1 * t34 * t21;
      res__4_11 = -1680 * theta5__dot * t49;
      t52 = theta3__dot ^ 2;
      res__5_4 = 89600 * t52 * t31;
      t55 = t19 * t52;
      t60 = t21 * t52;
      res__5_5 = t12 * (-960 * t9 - 89600 * t55 - 143520 * t16) - 143520 * t22 * (t9 + 0.560e3 / 0.897e3 * t60 - 0.2e1 / 0.299e3 * t16);
      t69 = t12 * t34 - 0.2e1 / 0.299e3 * t37 * t22;
      res__5_6 = 143520 * t8 * t69;
      res__5_9 = -179200 * theta3__dot * t44;
      t76 = t12 * t37 + 0.299e3 / 0.2e1 * t34 * t22;
      res__5_11 = -1920 * theta5__dot * t76;
      res__6_4 = 125580 * t52 * t40;
      res__6_5 = 143520 * t11 * t69;
      res__6_6 = t15 * (960 * t23 - 125580 * t55 + 840 * t60 - 143520 * t13) - 143520 * t7 * (t23 + 0.7e1 / 0.1196e4 * t55 + 0.7e1 / 0.8e1 * t60 + 0.2e1 / 0.299e3 * t13);
      res__6_9 = 1680 * theta3__dot * t49;
      res__6_10 = 1920 * theta4__dot * t76;
      t95 = cos(theta1);
      res__7_1 = t95 * t1;
      t96 = sin(theta1);
      t97 = theta1__dot * t96;
      t98 = s + L0;
      t99 = t98 * t97;
      t100 = s__dot * t95;
      res__7_2 = -theta1__dot * (t99 - 2 * t100);
      t104 = sin(theta2);
      res__7_3 = -L1 * t104 * t52;
      res__7_7 = 2 * t97;
      t107 = t95 * theta1__dot;
      t108 = t98 * t107;
      t109 = s__dot * t96;
      res__7_8 = 2 * t108 + 2 * t109;
      t111 = cos(theta2);
      res__7_9 = 2 * L1 * t111 * theta3__dot;
      res__8_1 = t96 * t1;
      res__8_2 = theta1__dot * (t108 + 2 * t109);
      res__8_3 = L1 * t111 * t52;
      res__8_7 = -2 * t107;
      res__8_8 = 2 * t99 - 2 * t100;
      res__8_9 = 2 * L1 * t104 * theta3__dot;
      res__10_4 = -L2 * t21 * t11;
      res__10_5 = -L4 * t22 * t8;
      t125 = theta2__dot ^ 2;
      res__10_6 = -L5 * t7 * t125;
      res__10_10 = 2 * L2 * t19 * theta4__dot;
      res__10_11 = 2 * L4 * t12 * theta5__dot;
      res__10_12 = 2 * L5 * t15 * theta2__dot;
      res__11_4 = L2 * t19 * t11;
      res__11_5 = L4 * t12 * t8;
      res__11_6 = L5 * t15 * t125;
      res__11_10 = 2 * L2 * t21 * theta4__dot;
      res__11_11 = 2 * L4 * t22 * theta5__dot;
      res__11_12 = 2 * L5 * t7 * theta2__dot;
      
      % store on output
      res__JbigRHS = zeros(11,12);
      res__JbigRHS(1,1) = res__1_1;
      res__JbigRHS(1,8) = res__1_8;
      res__JbigRHS(2,1) = res__2_1;
      res__JbigRHS(2,7) = res__2_7;
      res__JbigRHS(2,8) = res__2_8;
      res__JbigRHS(4,4) = res__4_4;
      res__JbigRHS(4,5) = res__4_5;
      res__JbigRHS(4,6) = res__4_6;
      res__JbigRHS(4,10) = res__4_10;
      res__JbigRHS(4,11) = res__4_11;
      res__JbigRHS(5,4) = res__5_4;
      res__JbigRHS(5,5) = res__5_5;
      res__JbigRHS(5,6) = res__5_6;
      res__JbigRHS(5,9) = res__5_9;
      res__JbigRHS(5,11) = res__5_11;
      res__JbigRHS(6,4) = res__6_4;
      res__JbigRHS(6,5) = res__6_5;
      res__JbigRHS(6,6) = res__6_6;
      res__JbigRHS(6,9) = res__6_9;
      res__JbigRHS(6,10) = res__6_10;
      res__JbigRHS(7,1) = res__7_1;
      res__JbigRHS(7,2) = res__7_2;
      res__JbigRHS(7,3) = res__7_3;
      res__JbigRHS(7,7) = res__7_7;
      res__JbigRHS(7,8) = res__7_8;
      res__JbigRHS(7,9) = res__7_9;
      res__JbigRHS(8,1) = res__8_1;
      res__JbigRHS(8,2) = res__8_2;
      res__JbigRHS(8,3) = res__8_3;
      res__JbigRHS(8,7) = res__8_7;
      res__JbigRHS(8,8) = res__8_8;
      res__JbigRHS(8,9) = res__8_9;
      res__JbigRHS(10,4) = res__10_4;
      res__JbigRHS(10,5) = res__10_5;
      res__JbigRHS(10,6) = res__10_6;
      res__JbigRHS(10,10) = res__10_10;
      res__JbigRHS(10,11) = res__10_11;
      res__JbigRHS(10,12) = res__10_12;
      res__JbigRHS(11,4) = res__11_4;
      res__JbigRHS(11,5) = res__11_5;
      res__JbigRHS(11,6) = res__11_6;
      res__JbigRHS(11,10) = res__11_10;
      res__JbigRHS(11,11) = res__11_11;
      res__JbigRHS(11,12) = res__11_12;

    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res__bigM = bigM( self, t, vars__ )

      g  = self.gravity;
      m1 = self.m1;
      m2 = self.m2;
      L1 = self.L1;
      L2 = self.L2;

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta2 = vars__(3);
      theta3 = vars__(4);
      theta4 = vars__(5);
      theta5 = vars__(6);
      s__dot = vars__(7);
      theta1__dot = vars__(8);
      theta3__dot = vars__(9);
      theta4__dot = vars__(10);
      theta5__dot = vars__(11);
      theta2__dot = vars__(12);

      % evaluate function
      res__1_1 = 5;
      res__1_2 = -6;
      res__1_7 = cos(theta1);
      res__1_8 = sin(theta1);
      res__2_1 = -6;
      t1 = s ^ 2;
      res__2_2 = 5 * t1 + 1602 * s + 0.256689e6 / 0.2e1;
      t4 = s + L0;
      res__2_7 = -t4 * res__1_8;
      res__2_8 = t4 * res__1_7;
      t6 = sin(theta2);
      res__3_7 = -L1 * t6;
      t8 = cos(theta2);
      res__3_8 = L1 * t8;
      res__3_9 = -1;
      res__4_3 = 98000;
      t9 = sin(theta3);
      t10 = sin(theta4);
      t12 = cos(theta3);
      t13 = cos(theta4);
      res__4_4 = 89600 * t10 * t9 + 89600 * t13 * t12;
      t16 = cos(theta5);
      t18 = sin(theta5);
      t23 = t16 - 0.299e3 / 0.2e1 * t18;
      res__4_5 = t12 * (125580 * t16 + 840 * t18) - 840 * t23 * t9;
      res__4_9 = 1;
      res__4_10 = -L2 * t9;
      res__4_11 = L2 * t12;
      res__5_3 = res__4_4;
      res__5_4 = 89603;
      res__5_5 = t13 * (960 * t18 + 143520 * t16) - 960 * t23 * t10;
      res__5_10 = -L4 * t10;
      res__5_11 = L4 * t13;
      res__6_3 = res__4_5;
      res__6_4 = res__5_5;
      res__6_5 = 268216;
      res__6_10 = -L5 * t18;
      res__6_11 = L5 * t16;
      res__7_1 = res__1_7;
      res__7_2 = res__2_7;
      res__7_3 = res__3_7;
      res__8_1 = res__1_8;
      res__8_2 = res__2_8;
      res__8_3 = res__3_8;
      res__9_3 = -1;
      res__9_4 = 1;
      res__10_4 = res__4_10;
      res__10_5 = res__5_10;
      res__10_6 = res__6_10;
      res__11_4 = res__4_11;
      res__11_5 = res__5_11;
      res__11_6 = res__6_11;
      
      % store on output
      res__bigM = zeros(11,11);
      res__bigM(1,1) = res__1_1;
      res__bigM(1,2) = res__1_2;
      res__bigM(1,7) = res__1_7;
      res__bigM(1,8) = res__1_8;
      res__bigM(2,1) = res__2_1;
      res__bigM(2,2) = res__2_2;
      res__bigM(2,7) = res__2_7;
      res__bigM(2,8) = res__2_8;
      res__bigM(3,7) = res__3_7;
      res__bigM(3,8) = res__3_8;
      res__bigM(3,9) = res__3_9;
      res__bigM(4,3) = res__4_3;
      res__bigM(4,4) = res__4_4;
      res__bigM(4,5) = res__4_5;
      res__bigM(4,9) = res__4_9;
      res__bigM(4,10) = res__4_10;
      res__bigM(4,11) = res__4_11;
      res__bigM(5,3) = res__5_3;
      res__bigM(5,4) = res__5_4;
      res__bigM(5,5) = res__5_5;
      res__bigM(5,10) = res__5_10;
      res__bigM(5,11) = res__5_11;
      res__bigM(6,3) = res__6_3;
      res__bigM(6,4) = res__6_4;
      res__bigM(6,5) = res__6_5;
      res__bigM(6,10) = res__6_10;
      res__bigM(6,11) = res__6_11;
      res__bigM(7,1) = res__7_1;
      res__bigM(7,2) = res__7_2;
      res__bigM(7,3) = res__7_3;
      res__bigM(8,1) = res__8_1;
      res__bigM(8,2) = res__8_2;
      res__bigM(8,3) = res__8_3;
      res__bigM(9,3) = res__9_3;
      res__bigM(9,4) = res__9_4;
      res__bigM(10,4) = res__10_4;
      res__bigM(10,5) = res__10_5;
      res__bigM(10,6) = res__10_6;
      res__bigM(11,4) = res__11_4;
      res__bigM(11,5) = res__11_5;
      res__bigM(11,6) = res__11_6;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res__JbigETA = JbigETA( self, t, vars__, mu )

      mu1 = mu(1);
      mu2 = mu(2);
      mu3 = mu(3);
      mu4 = mu(4);
      mu5 = mu(5);
      mu6 = mu(6);
      mu7 = mu(7);
      mu8 = mu(8);

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta2 = vars__(3);
      theta3 = vars__(4);
      theta4 = vars__(5);
      theta5 = vars__(6);
      s__dot = vars__(7);
      theta1__dot = vars__(8);
      theta3__dot = vars__(9);
      theta4__dot = vars__(10);
      theta5__dot = vars__(11);
      &theta;2__dot = vars__(12);

      % evaluate function
      t1 = sin(theta1);
      t2 = mu7 * t1;
      t3 = cos(theta1);
      t4 = mu8 * t3;
      res__1_2 = -t2 + t4;
      res__2_1 = (10 * s + 1602) * mu2 - t2 + t4;
      t8 = -s - L0;
      res__2_2 = (mu7 * t3 + mu8 * t1) * t8;
      t12 = sin(theta2);
      t14 = cos(theta2);
      res__3_3 = -(mu7 * t14 + mu8 * t12) * L1;
      t19 = cos(theta5);
      t20 = mu5 * t19;
      t22 = sin(theta4);
      t23 = t22 * mu4;
      t25 = sin(theta5);
      t26 = mu5 * t25;
      t29 = cos(theta3);
      t31 = sin(theta3);
      t33 = cos(theta4);
      t34 = mu4 * t33;
      res__4_4 = t29 * (-L2 * mu10 - 840 * t20 + 89600 * t23 + 125580 * t26) - (L2 * mu11 + 125580 * t20 + 840 * t26 + 89600 * t34) * t31;
      t42 = t22 * t29 - t31 * t33;
      res__4_5 = -89600 * mu4 * t42;
      t46 = t19 - 0.299e3 / 0.2e1 * t25;
      t49 = t19 + 0.2e1 / 0.299e3 * t25;
      t52 = t29 * t46 + 0.299e3 / 0.2e1 * t49 * t31;
      res__4_6 = 840 * mu5 * t52;
      res__5_4 = 89600 * mu3 * t42;
      t57 = mu3 * t31;
      t63 = mu3 * t29;
      res__5_5 = t33 * (-L4 * mu10 - 960 * t20 + 143520 * t26 + 89600 * t57) - (L4 * mu11 + 143520 * t20 + 960 * t26 + 89600 * t63) * t22;
      t72 = t33 * t46 + 0.299e3 / 0.2e1 * t49 * t22;
      res__5_6 = 960 * mu5 * t72;
      res__6_4 = -840 * mu3 * t52;
      res__6_5 = -960 * mu4 * t72;
      res__6_6 = t19 * (-L5 * mu10 + 143520 * t23 + 960 * t34 + 125580 * t57 + 840 * t63) - (L5 * mu11 - 960 * t23 + 143520 * t34 - 840 * t57 + 125580 * t63) * t25;
      res__7_1 = -mu2 * t1;
      t94 = -t8;
      res__7_2 = -mu2 * t94 * t3 - mu1 * t1;
      res__7_3 = -mu3 * L1 * t14;
      res__8_1 = mu2 * t3;
      res__8_2 = -mu2 * t94 * t1 + mu1 * t3;
      res__8_3 = -mu3 * L1 * t12;
      res__10_4 = -mu4 * L2 * t29;
      res__10_5 = -mu5 * L4 * t33;
      res__10_6 = -mu6 * L5 * t19;
      res__11_4 = -mu4 * L2 * t31;
      res__11_5 = -mu5 * L4 * t22;
      res__11_6 = -mu6 * L5 * t25;
      
      % store on output
      res__JbigETA = zeros(11,12);
      res__JbigETA(1,2) = res__1_2;
      res__JbigETA(2,1) = res__2_1;
      res__JbigETA(2,2) = res__2_2;
      res__JbigETA(3,3) = res__3_3;
      res__JbigETA(4,4) = res__4_4;
      res__JbigETA(4,5) = res__4_5;
      res__JbigETA(4,6) = res__4_6;
      res__JbigETA(5,4) = res__5_4;
      res__JbigETA(5,5) = res__5_5;
      res__JbigETA(5,6) = res__5_6;
      res__JbigETA(6,4) = res__6_4;
      res__JbigETA(6,5) = res__6_5;
      res__JbigETA(6,6) = res__6_6;
      res__JbigETA(7,1) = res__7_1;
      res__JbigETA(7,2) = res__7_2;
      res__JbigETA(7,3) = res__7_3;
      res__JbigETA(8,1) = res__8_1;
      res__JbigETA(8,2) = res__8_2;
      res__JbigETA(8,3) = res__8_3;
      res__JbigETA(10,4) = res__10_4;
      res__JbigETA(10,5) = res__10_5;
      res__JbigETA(10,6) = res__10_6;
      res__JbigETA(11,4) = res__11_4;
      res__JbigETA(11,5) = res__11_5;
      res__JbigETA(11,6) = res__11_6;

    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = DoublePendulum8EQ( data, eta, omega )
      % call the constructor of the basic class
      neq  = 8;
      ninv = 4;
      self@DAC_ODEclass('DoublePendulum8EQ',neq,ninv);
      % setup of the parmater of the ODE

      self.gravity = data.gravity;
      self.L1      = data.L1;
      self.L2      = data.L2;
      self.m1      = data.m1;
      self.m2      = data.m2;
      self.eta     = eta;
      self.omega   = omega;
      self.npos    = 4;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setBaumgarte( self, eta, omega )
      self.eta   = eta;
      self.omega = omega;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res__rhs = f( self, t, vars__ )
      npos = self.npos;
      %
      %       / f(q,v,t) \
      % RHS = |          |
      %       \ g(q,v,t) /
      %
      RHS = self.bigRHS( t, vars__ );
      %
      %
      %        / M(q,v,t)  Phi_q(q,t)^T \
      % BIGM = |                        |
      %        \ Phi_q(q,t)     0       /
      %
      BIGM = self.bigM( t, vars__ );
      SOL  = BIGM\RHS;
    
      % evaluate function
      res__rhs = zeros(2*npos,1);
      res__rhs(1:npos)        = vars__(npos+1:2*npos);
      res__rhs(npos+1:2*npos) = SOL(1:npos);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function jac = DfDx( self, t, vars__ )
      %
      %       / f(q,v,t) \
      % RHS = |          |
      %       \ g(q,v,t) /
      %
      RHS = self.bigRHS( t, vars__ );
      %
      %
      %        / M(q,v,t)  Phi_q(q,t)^T \
      % BIGM = |                        |
      %        \ Phi_q(q,t)     0       /
      %
      BIGM = self.bigM( t, vars__ );

      SOL  = BIGM\RHS;
      JETA = self.JbigETA( t, vars__, SOL );
      JRHS = self.JbigRHS( t, vars__ );
      JM   = BIGM\(JRHS-JETA);
      
      % combine all jacobians
      npos = self.npos;
      jac  = [ zeros(npos,npos) eye(npos); JM(1:npos,:) ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res__h = h( self, t, vars__ )
      g  = self.gravity;
      m1 = self.m1;
      m2 = self.m2;
      L1 = self.L1;
      L2 = self.L2;

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta2 = vars__(3);
      theta3 = vars__(4);
      theta4 = vars__(5);
      theta5 = vars__(6);
      s__dot = vars__(7);
      theta1__dot = vars__(8);
      theta3__dot = vars__(9);
      theta4__dot = vars__(10);
      theta5__dot = vars__(11);
      theta2__dot = vars__(12);
      s__dot__d = vars__(13);
      theta1__dot__d = vars__(14);
      theta3__dot__d = vars__(15);
      theta4__dot__d = vars__(16);
      theta5__dot__d = vars__(17);
      theta2__dot__d = vars__(18);

      % evaluate function
      res__1 = s;
      res__2 = theta1;
      res__3 = theta2;
      res__4 = theta3;
      res__5 = theta4;
      res__6 = theta5;
      res__7 = lambda1(t);
      res__8 = lambda2(t);
      res__9 = lambda3(t);
      res__10 = lambda4(t);
      res__11 = lambda5(t);
      res__12 = s__dot;
      res__13 = theta1__dot;
      res__14 = theta3__dot;
      res__15 = theta4__dot;
      res__16 = theta5__dot;
      res__17 = theta2__dot;
      

      % store on output
      res__h = zeros(17,1);
      res__h(1) = res__1;
      res__h(2) = res__2;
      res__h(3) = res__3;
      res__h(4) = res__4;
      res__h(5) = res__5;
      res__h(6) = res__6;
      res__h(7) = res__7;
      res__h(8) = res__8;
      res__h(9) = res__9;
      res__h(10) = res__10;
      res__h(11) = res__11;
      res__h(12) = res__12;
      res__h(13) = res__13;
      res__h(14) = res__14;
      res__h(15) = res__15;
      res__h(16) = res__16;
      res__h(17) = res__17;

    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res__DhDx = DhDx( self, t, vars__ )
      g  = self.gravity;
      m1 = self.m1;
      m2 = self.m2;
      L1 = self.L1;
      L2 = self.L2;

      % extract states
      s = vars__(1);
      theta1 = vars__(2);
      theta2 = vars__(3);
      theta3 = vars__(4);
      theta4 = vars__(5);
      theta5 = vars__(6);
      s__dot = vars__(7);
      theta1__dot = vars__(8);
      theta3__dot = vars__(9);
      theta4__dot = vars__(10);
      theta5__dot = vars__(11);
      theta2__dot = vars__(12);
      s__dot__d = vars__(13);
      theta1__dot__d = vars__(14);
      theta3__dot__d = vars__(15);
      theta4__dot__d = vars__(16);
      theta5__dot__d = vars__(17);
      theta2__dot__d = vars__(18);

      % evaluate function
      t1 = sin(theta1);
      t2 = theta1__dot * t1;
      res__1_1 = -t2;
      t3 = cos(theta1);
      t5 = s + L0;
      res__1_2 = -t5 * t3 * theta1__dot - s__dot * t1;
      t8 = cos(theta2);
      res__1_3 = -L1 * t8 * theta3__dot;
      res__1_7 = t3;
      res__1_8 = -t5 * t1;
      t12 = sin(theta2);
      res__1_9 = -L1 * t12;
      res__2_1 = res__1_7 * theta1__dot;
      res__2_2 = s__dot * res__1_7 - t5 * t2;
      res__2_3 = -L1 * t12 * theta3__dot;
      res__2_7 = t1;
      res__2_8 = t5 * res__1_7;
      res__2_9 = L1 * t8;
      res__3_9 = -1;
      res__3_10 = 1;
      t18 = cos(theta3);
      res__4_4 = -L2 * t18 * theta4__dot;
      t21 = cos(theta4);
      res__4_5 = -L4 * t21 * theta5__dot;
      t24 = cos(theta5);
      res__4_6 = -L5 * t24 * theta2__dot;
      t27 = sin(theta3);
      res__4_10 = -L2 * t27;
      t29 = sin(theta4);
      res__4_11 = -L4 * t29;
      t31 = sin(theta5);
      res__4_12 = -L5 * t31;
      res__5_4 = -L2 * t27 * theta4__dot;
      res__5_5 = -L4 * t29 * theta5__dot;
      res__5_6 = -L5 * t31 * theta2__dot;
      res__5_10 = L2 * t18;
      res__5_11 = L4 * t21;
      res__5_12 = L5 * t24;
      
      % store on output
      res__DhDx = zeros(5,18);
      res__DhDx(1,1) = res__1_1;
      res__DhDx(1,2) = res__1_2;
      res__DhDx(1,3) = res__1_3;
      res__DhDx(1,7) = res__1_7;
      res__DhDx(1,8) = res__1_8;
      res__DhDx(1,9) = res__1_9;
      res__DhDx(2,1) = res__2_1;
      res__DhDx(2,2) = res__2_2;
      res__DhDx(2,3) = res__2_3;
      res__DhDx(2,7) = res__2_7;
      res__DhDx(2,8) = res__2_8;
      res__DhDx(2,9) = res__2_9;
      res__DhDx(3,9) = res__3_9;
      res__DhDx(3,10) = res__3_10;
      res__DhDx(4,4) = res__4_4;
      res__DhDx(4,5) = res__4_5;
      res__DhDx(4,6) = res__4_6;
      res__DhDx(4,10) = res__4_10;
      res__DhDx(4,11) = res__4_11;
      res__DhDx(4,12) = res__4_12;
      res__DhDx(5,4) = res__5_4;
      res__DhDx(5,5) = res__5_5;
      res__DhDx(5,6) = res__5_6;
      res__DhDx(5,10) = res__5_10;
      res__DhDx(5,11) = res__5_11;
      res__DhDx(5,12) = res__5_12;
    end

    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, t, Z )
      DoublePendulumPlot( t, Z(1), Z(2), Z(3), Z(4), self.L1, self.L2 );
    end
  end
end
