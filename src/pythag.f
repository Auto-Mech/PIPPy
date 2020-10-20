      double precision function pythag(a,b)
      double precision a,b
      double precision p,q,r,s,t
      p = dmax1(abs(a),abs(b))
      q = dmin1(abs(a),abs(b))
      if (q .eq. 0.0e0) go to 20
   10 continue
         r = (q/p)**2
         t = 4.0e0 + r
         if (t .eq. 4.0e0) go to 20
         s = r/t
         p = p + 2.0e0*p*s
         q = q*s
      go to 10
   20 pythag = p
      return
      end

