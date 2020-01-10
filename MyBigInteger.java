import java.util.ArrayList;

public class MyBigInteger  implements Comparable<MyBigInteger> {   //  We will cover what "implements Comparable" means later

	/*
	 *   
	 *  If the number has N digits,  then it is represented is a polynomial:
	 *  coefficients[N-1] base^{N-1} + ...  coefficients[1] base^{1} +  coefficients[0] 
	 *  where base has a particular value and the coefficients are in {0, 1, ...  base - 1}
	 *  
	 *  For any base and any positive integer, the representation of that positive 
	 *  integer as a sum of powers of that base is unique.  
	 *  
	 *  We require that the coefficient of the largest power is non-zero.  
	 *  For example,  '354' is a valid representation (which we call "three hundred fifty four") 
	 *  but '0354' is not.  
	 * 
	 */
	
	private int	base;       

	private ArrayList<Integer>  coefficients;

	//  Constructors

	MyBigInteger(int base){
		
		//  If no string argument is given, then this constructor constructs an empty list of coefficients.
		
		this.base = base;
		coefficients = new ArrayList<Integer>();
	}

	/*       
	 *   The constructor builds a ArrayList of Integer objects where the integers are in [0,base), i.e.  0, 1, ..., base - 1.
	 *   Note that the argument numberAsString only represents a base 10 number when the base is 10.
	 */
	
	MyBigInteger(String numberAsString,  int base) throws Exception{
		int i;
		this.base = base;
		coefficients = new ArrayList<Integer>();
		if ((base < 2) || (base > 10)){
			System.out.println("constructor error:  base must be between 2 and 10 (inclusive)");
			throw new Exception();
		}

		int len = numberAsString.length();
		for (int indx = 0; indx < len; indx++){  
			i = numberAsString.charAt(indx);
			
			/*
			 *  ascii value of symbol '0' is 48, symbol '1' is 49, etc.
			 *  So,  to get the numerical value of '2',  we subtract 
             *  the character value of '0' (48) from the character value of '2' (50).
			 */
            
			if ( (i >= 48) && (i - 48 < base))   				
				coefficients.add(0,  Integer.valueOf(i-48) );
			else{
				System.out.println("constructor error:  all coefficients should be non-negative and less than base");
				throw new Exception();
			}			
		}
	}
	
	/*
	 *   This Constructor allows the input to be an int.
	 *   The maximum value can therefore be 2^31-1.     
	 *   Note that the method assumes the input i is written in base 10.
	 *   It then converts that value to a particular base.
	 */
	
	MyBigInteger(int i,  int base) throws Exception{
		this.base = base;
		coefficients = new ArrayList<Integer>();
		
		if (i < 0) {
			System.out.println("constructor error: we are only dealing with non-negative integers");
			throw new Exception();
		}
		else if (i == 0)
			coefficients.add( Integer.valueOf(0) );
		else {
			int m = i;
			while (m > 0) {
				coefficients.add( Integer.valueOf(m % base) );
				m /= base;
			}
		}
	}
	
	public int getBase()
	{
		return base;
	}
	
	/*  
	 * Normally we would not need this method, since coefficients are meant to be hidden from the user.
	 * The reason for having the method is that the Grader can use it to examine the number.
	 * Another way to do this would have been to make the field public. 
	 */
	
	public ArrayList<Integer>  getCoefficients()
	{
		return coefficients;
	}

	/*
	 *   The plus method computes this.plus(b), that is, a+b where 'this' is a.
	 */

	public MyBigInteger plus( MyBigInteger  second) throws Exception{
				
		//  initialize the result as an empty list of coefficients
		
		MyBigInteger result = new MyBigInteger( this.base );

		if (this.base != second.base){
			throw new Exception("Bases must be the same in an addition");
		}

		/*   
		 * The plus method must not affect the numbers themselves. 
		 * So let's just work  with a copy (a clone) of the numbers. 
		 */

		MyBigInteger  firstClone  = this.clone();
		MyBigInteger  secondClone = second.clone();

		/*
		 *   If the two numbers have a different polynomial order    
		 *   then pad the smaller one with zero coefficients.
		 */
		
		int   diff = firstClone.coefficients.size() - second.coefficients.size();
		while (diff < 0){  // second is bigger                                                     		
			firstClone.coefficients.add(0);			        
			diff++;
		}
		while (diff > 0){  //  this is bigger
			secondClone.coefficients.add(0);
			diff--;
		}

		/*
		 *   'firstClone' and 'secondClone' have the same size.  We add the coefficients
		 *   term by term.    If the last coefficient yields a carry, then we add 
         *   one more term with the carry.	
         */
		
		int tmp;
		int carry = 0;
		
		for (int i=0; i < firstClone.coefficients.size(); i++  ){
			tmp = firstClone.coefficients.get(i) + secondClone.coefficients.get(i) + carry;
			result.coefficients.add( tmp % base ) ;
			carry = tmp / base; 
		}
		if (carry > 0)
			result.coefficients.add(carry);    //   carry always would have the value 1 in this case.
		
		return result;		
	}

	/*   
	 *    Slow multiplication algorithm.
	 *    'this' is the multiplicand i.e. a*b = a+a+a+...+a (b times) where a is multiplicand and b is multiplier
	 */
	
	public MyBigInteger slowTimes( MyBigInteger  multiplier) throws Exception{
		
		if (this.base != multiplier.base){
			throw new Exception("Bases must be the same in slowTimes");
		}

		MyBigInteger prod  = new MyBigInteger(0, this.base);
		MyBigInteger one   = new MyBigInteger(1, this.base);  // This is used to increment the counter.
		
		for (MyBigInteger counter = new MyBigInteger(0, this.base);  counter.compareTo(multiplier) < 0;  counter = counter.plus(one) ){
			prod = prod.plus(this);
		}
		return prod;
	}
	
	/*
	 *    The times method is NOT be the same as what you learned in grade school which 
	 *    uses a temporary 2D table with space proportional to the square of 
	 *    the number of coefficients in the operands, i.e.  N^2.   The method here accumulates 
	 *    the result by adding each row of the table as it is computed. 
	 *    This method uses space that is proportional to the number of coefficients N.    
	 *    The multiplication algorithm will still take time O(N^2) however.   
	 *    
	 *    The  method computes this.times(b) where 'this' is a.
	 */
	
	public MyBigInteger times( MyBigInteger multiplicand) throws Exception{
		
		if (this.base != multiplicand.base){
			throw new Exception("Bases must be the same in times()");
		}
		
		//  anything times 0 is 0,  so if either 'this' or multiplicand are 0 then return 0.
		MyBigInteger zero   = new MyBigInteger(0,this.base);   
		if ((this.compareTo(zero) == 0) || (multiplicand.compareTo(zero) == 0))
			return zero;
		
		//  initialize product as an empty list of coefficients
		
		MyBigInteger product = new MyBigInteger( this.base );

		/*
		 *           multiplicand
		 *          x  multiplier  (this)
		 *        ---------------
		 *            
		 *   Note we use helper methods.
		 */

		MyBigInteger row;
		for (int i=0; i < this.coefficients.size(); i++){
			row = multiplicand.timesSingleDigit( this.coefficients.get(i) ).timesBaseToThePower(i);
			product = product.plus( row  );
		}
		return product;
	}

	/*
	 *    'this' (the caller) will be the multiplicand.   
	 */
	
	private MyBigInteger timesSingleDigit( int  singleDigit) throws Exception{

		//  Assumes that 0 <= singleDigit < base. 
		
		if ((singleDigit >= this.base) || (singleDigit < 0)){
			throw new Exception("Single digit must be in {0, ..., this.base - 1}  in timeSingleDigit");
		}
		
		//  if we multiply 'this' by 0, then the result is 0
		if (singleDigit == 0) {
			return new MyBigInteger(0, this.base);   
		}
		
		int carry = 0;
		int tmp;

		//  initialize prod as an empty list of coefficients
	
		MyBigInteger prod = new MyBigInteger(this.base);

		//  multiply the single digit by each of the digits in this.coefficients
		//  This is like a row in the grade school calculation.
		
		for (int i=0; i < this.coefficients.size(); i++  ){
			tmp = this.coefficients.get(i) * singleDigit + carry;
			prod.coefficients.add( tmp % this.base ) ;
			carry = tmp / this.base; 
		}
		if (carry > 0)
			prod.coefficients.add(carry);
		return prod;
	}
	
	/*
	 *   The minus method computes this.minus(b) where 'this' is a, and a > b.
	 *   If a < b, then it throws an exception.
	 *   
	 *   The solution below uses variable names  first and second rather than a and b.
	 *   So we are computing first.minus(second)
	 */
	
	public MyBigInteger  minus(MyBigInteger second) throws Exception{

		//  initialize the result as an empty list of coefficients
		
		MyBigInteger  result = new MyBigInteger(this.base);

		if (this.base != second.base){
			throw new Exception("Bases must be the same in minus method.");
		}
		/*
		 *    The minus method is not supposed to change the numbers. 
		 *    But the grade school algorithm sometimes requires us to "borrow"
		 *    from a higher coefficient to a lower one.   So we work
		 *    with a copy (a clone) instead.
		 */

		MyBigInteger  first = this.clone();

		//   Verify 'this' >= second. 
		 
		if (this.compareTo(second) < 0){
			throw new Exception("a.minus(b) requires that a >= b");
		}

		int i = 0;       //  coefficient position 
		int diffCoef;    //  compute the difference of two coefficients   this[i] - second[i]
	
		while (i < first.coefficients.size()){
		
			/*
			 *   Start from the least significant digit.
			 *   For each i, check if 'second' has a term at position i. 
			 *   If yes, then need to take the difference at this position.    That's where the bulk of work is.
			 *   If no, then the coefficient at second is treated as 0 and just copy the i'th coefficient from first to the result. 
			 */
			
			if (i < second.coefficients.size()){
				
				diffCoef = first.coefficients.get(i) - second.coefficients.get(i);		

				if (diffCoef >= 0)
					result.coefficients.add(Integer.valueOf(diffCoef));
				else {
					/*
					 *   .. then we need to borrow from the next coefficient.
					 *   But if the next coefficient holds a '0', then we need a sequence of borrows.
					 *   
					 *   For example, suppose we are computing 30001 - 6.   Then we essentially treat
					 *   30001 as  11 + 90 + 900 + 9000 + 2000.     
					 *   So we would want to rewrite the list of coefficients as (11, 9, 9, 9, 2).
					 */
					int j = i;
					while (first.coefficients.get(j+1) == 0){
						first.coefficients.set(j+1, base-1);
						j++;
					}
					first.coefficients.set(j+1, Integer.valueOf( first.coefficients.get(j+1) - 1 ));

					//  Finally, compute the result for position i
					result.coefficients.add(Integer.valueOf(base + diffCoef));
				}
			}
			else
				result.coefficients.add(Integer.valueOf(first.coefficients.get(i)));
			i++;
		}

		/*  
		 *  In the case of say  100-98, we will end up with 002.  
		 *  Remove all the leading 0's of the result.
		 */
		
		while ((result.coefficients.size() > 1) & 
				(result.coefficients.get(result.coefficients.size()-1) == 0)){
			result.coefficients.remove(result.coefficients.size()-1);
		}
		return result;	
	}

	/*   
	 *    Slow division algorithm is repeated subtraction, mentioned in lecture 1.
	 */

	public MyBigInteger slowdividedBy( MyBigInteger  divisor) throws Exception{

		if (this.base != divisor.base){
			throw new Exception("Bases must be the same in slowdividedBy method");
		}
		
		MyBigInteger one = new MyBigInteger(1,base);
		MyBigInteger quotient = new MyBigInteger(0,base);
		MyBigInteger remainder = this.clone();
		while ( remainder.compareTo(divisor) >= 0 ){
			remainder = remainder.minus(divisor);
			quotient = quotient.plus(one);
		}
		return quotient;
	}
	
	/*  
	 *  The dividedBy method divides 'this' by 'divisor' i.e. this.dividedBy(divisor)
	 *   It returns the quotient and ignores the remainder.
	 */
	
	public MyBigInteger dividedBy( MyBigInteger  divisor ) throws Exception{
		
		//   ADD YOUR CODE BELOW HERE
		
		if (this.base != divisor.base){
			throw new Exception("Bases must be the same in dividedBy method");
		}
		
		MyBigInteger quotient = new MyBigInteger(this.base); // your answer with the same base, empty list of coefficients
		MyBigInteger dividend = this.clone(); // the integer that will be divided by the divisor
		MyBigInteger tmp = new MyBigInteger(this.base); // empty list of coefficients!
		MyBigInteger zeroInteger = new MyBigInteger (0, this.base); // just the zero integer
		MyBigInteger oneInteger = new MyBigInteger(1, this.base); // just the one integer
		
		if (this.compareTo(divisor) == 0) { // if the dividend and the divisor are equal,
			 return oneInteger; // the answer is simply 1.
		} else if (this.compareTo(divisor) < 0) { // if the dividend is smaller than the divisor,
			return zeroInteger; // the answer is simply 0.
		}
		
		int differenceInSize;
		differenceInSize = (dividend.coefficients.size()) - (divisor.coefficients.size());
		//System.out.print("difference in size is: "); System.out.println(differenceInSize); this is for debugging
		
		if (differenceInSize > 1) {
			for (int i = differenceInSize - 1; i < dividend.coefficients.size(); i++) {
				tmp.coefficients.add(dividend.coefficients.get(i));
				//System.out.print("tmp coefficients:"); System.out.println(tmp.getCoefficients());
				}
			for (int a = dividend.coefficients.size() - 1; a >= differenceInSize - 1; a--) {
				dividend.coefficients.remove(a);
				//System.out.print("dividend coefficients: "); System.out.println(dividend.getCoefficients());
			}
		} else {
			for (int k = differenceInSize; k < dividend.coefficients.size(); k++) {
				tmp.coefficients.add(dividend.coefficients.get(k));
				//System.out.print("tmp coefficients:"); System.out.println(tmp.getCoefficients());
			}
			for (int b = dividend.coefficients.size() - 1; b >= differenceInSize; b--) {
				dividend.coefficients.remove(b);
			}
		}
				
        while (dividend.coefficients.size() >= 0) {
            MyBigInteger slowdivision = tmp.slowdividedBy(divisor);
            for (int x = slowdivision.coefficients.size() - 1; x >= 0; x--) {
            	int h = 0;
            	quotient.coefficients.add(h, slowdivision.coefficients.get(x));
            	h++;
            }
            
            MyBigInteger multiplication = slowdivision.times(divisor);
            
            int firstvalue = multiplication.coefficients.get(multiplication.coefficients.size() - 1);
            
            if (firstvalue == 0) { // this is to remove leading 0s
            	while(firstvalue == 0) {
            		if (multiplication.coefficients.size() > 1) {
            			multiplication.coefficients.remove(firstvalue);
            		} else {
            			break;
            		}
            	}
            }
            
            int first = tmp.coefficients.get(tmp.coefficients.size() - 1);
            
            if (first == 0) { // this is to remove leading 0s
            	while(first == 0) {
            		if (tmp.coefficients.size() > 1) {
            			tmp.coefficients.remove(first);
            		} else {
            			break;
            		}
            	}
            }
            
            tmp = tmp.minus(multiplication);
            
            if (tmp.compareTo(zeroInteger) == 0) {
            	for (int r = tmp.coefficients.size() - 1; r >= 0; r--) {
            		tmp.coefficients.remove(r);
            	}
            }
            
            if (dividend.coefficients.size() == 0) {
                break;
            }
            
            tmp.coefficients.add(0, dividend.coefficients.remove(dividend.coefficients.size() - 1));
            
        }
				
		//  ADD YOUR CODE ABOVE HERE.

		return quotient;
	}
	
	/*  The convert method converts between two bases that are each in {2, ..., 10}
	 *  The convert method specifies the new base that you are converting into.
	 *  The subtlety of the convert method is that it uses operations in this.base
	 *  This can be tricky since if you are converting into a higher base then you 
	 *  have to consider that this.dividedBy() is performed in this.base. 
	 *  
	 */
	
	public MyBigInteger convert(int newBase) throws Exception {
		
		MyBigInteger remainder = this.clone();   
		MyBigInteger result = new MyBigInteger(newBase);
		
    	//   ADD YOUR CODE BELOW HERE
		
		// if current base equals to desired base, we return the current MyBigInteger
		
		if (this.base == newBase) {
			return this;
		}
		
		// if the current base does not equal the desired base, first, convert "this" to base 10 (decimal)
		
		MyBigInteger b = new MyBigInteger(this.base, 10);
		MyBigInteger s = new MyBigInteger(10);
		
		for (int i = remainder.coefficients.size() - 1; i >= 0; i--) {
			s = s.times(b);
			s = s.plus(new MyBigInteger(remainder.coefficients.get(i), 10));
		}
		
		// now, convert from base 10 to desired newBase
				
		MyBigInteger basebaby = new MyBigInteger (newBase, 10);
		while (s.coefficients.size() > 0) {
			MyBigInteger remainder2 = s.clone().mod(basebaby);
			result.coefficients.add(remainder2.coefficients.get(0));
			s = s.minus(remainder2);
			s = s.dividedBy(basebaby);
			if (s.coefficients.size() == 1 && s.coefficients.get(0) == 0) {
				s.coefficients.remove(0);
			}
		}
		
		//   ADD YOUR CODE ABOVE HERE
		
		return result;
	}
	
	public ArrayList<MyBigInteger>  primeFactors() throws Exception {
		
		ArrayList<MyBigInteger> factors = new ArrayList<MyBigInteger>();

		//  ADD YOUR CODE BELOW HERE
		
		MyBigInteger m = this.clone(); // set m to take the value of "this" MyBigInteger, m is a clone
		MyBigInteger twoInteger = new MyBigInteger(2, this.base);
		MyBigInteger threeInteger = new MyBigInteger(3, this.base);
		MyBigInteger zeroInteger = new MyBigInteger(0, this.base);
		
		// first, for efficiency purposes, add to the factors list the number of 2s that divide m	
		
		while (m.mod(twoInteger).compareTo(zeroInteger) == 0) { 
			factors.add(twoInteger);
			m = m.dividedBy(twoInteger);
		}
		
		// at this point, we know that m must be odd. Thus, we can skip one element (i = i + 2)
		
		for (MyBigInteger i = threeInteger; i.times(i).compareTo(m) <= 0 ; i = i.plus(twoInteger)) {
			while (m.mod(i).compareTo(zeroInteger) == 0) { // if the remainder of m mod i is equal to zero, then i is a prime factor of m 
				factors.add(i); // add MyBigInteger i (prime factor) to the arraylist "factors"
				m = m.dividedBy(i); // update the value of m by dividing m by i
			}
		}
		
		// If m is a prime number and m is greater than 2, do the following
		
		if (m.compareTo(twoInteger) == 1) {
			factors.add(m);
		}

		//  ADD YOUR CODE ABOVE HERE

		return factors;
	}

	//   ----------------   HELPER METHODS --------------------------
	
	/*  
	 *  The mod method divides 'this' by 'divisor' and returns the remainder.
	*/
	 
	public MyBigInteger mod( MyBigInteger  divisor ) throws Exception{
		
		if (this.base != divisor.base){
			throw new Exception("Bases must be the same in mod method");
		}
		return this.minus(this.dividedBy(divisor).times(divisor));		
	}		

	@Override
	public MyBigInteger  clone(){

		//  For technical reasons that don't interest us here, this method 
		//  has to be declared public (not private).

		MyBigInteger copy = new MyBigInteger(this.base);
		for (int i=0; i < this.coefficients.size(); i++){
			copy.coefficients.add( Integer.valueOf( this.coefficients.get(i) ) );
		}
		return copy;
	}
	
	/*
	 *  The a.compareTo(b) method returns -1 if a < b,  it returns 0 if a == b,  
	 *  and it returns 1 if a > b.
	 *  
 	 *  Assumes that numbers have valid representation  e.g.  no leading 0's.
	 */
	
	public int compareTo(MyBigInteger second) {

		//   if  this < second,  return -1  
		//   if  this > second,  return  1  
		//   otherwise they are equal and return 0
		
		if (this.base != second.base){
			throw new RuntimeException("Bases must be the same in compareTo method");
		}
		
		/*   
		 * Assume maximum degree coefficient is non-zero (except for case that one of the numbers is 0).   
		 * Then,  if two numbers have different maximum degree, it is easy to decide which is larger.
		 */

		int diff = this.coefficients.size() - second.coefficients.size();
		if (diff < 0)
			return -1;
		else if (diff > 0)
			return 1;
		else { 
			
			/*   
			 * If two numbers have the same maximum degree,  then it is a bit trickier
			 * to decide which number is larger.   You need to compare the coefficients,
			 * starting from the largest and working toward the smallest until you find
			 * coefficients that are not equal.
			 */
			
			boolean done = false;
			int i = this.coefficients.size() - 1;
			while (i >=0 && !done){
				diff = this.coefficients.get(i) - second.coefficients.get(i); 
				if (diff < 0){
					return -1;
				}
				else if (diff > 0)
					return 1;
				else{ 
					i--;
				}
			}
			return 0;    //   if all coefficients are the same,  so numbers are equal.
		}
	}

	/*  
	 *    computes  'this' * base^n  
	 *    Note that it modifies the number, so if you don't want that, then clone the number first.
	 */
	
	private MyBigInteger timesBaseToThePower(int n) throws Exception{
		
		if (n < 0){
			throw new Exception("timesBaseToThePower requires that the power (exponent) >= 0");
		}
		
		if (this.compareTo(new MyBigInteger(0,base)) == 0)
			return this;
		
		MyBigInteger thisShifted = new MyBigInteger(base);
		//  First make a list of n 0's.
		for (int i=0; i< n; i++){
			thisShifted.coefficients.add(Integer.valueOf(0));
		}
		
		//  Then add all the coefficents of 'this' after the n 0's.
		for (int i=0; i< this.coefficients.size(); i++){
			thisShifted.coefficients.add(this.coefficients.get(i));
		}
		return thisShifted;
	}

	/*
	 * 	Returns a string with coefficients in the reverse order which is the natural format for people to reading numbers,
	 *  i.e. people want to read  a[N-1], ... a[2] a[1] a[0]. 
	 */
	
	@Override
	public String toString(){	
		String s = new String(); 
		for (Integer coef : coefficients)     //  Java enhanced for loop
			s = coef.toString() + s ;        //   Append each successive coefficient.
		return "(" + s + ")_" + base;		
	}


	public static ArrayList<MyBigInteger> primesToN(int n, int base) throws Exception {
		
		//  Sieve of Eratosthenes algorithm
		// 
		//  Code here is modified from
		//  https://introcs.cs.princeton.edu/java/14array/PrimeSieve.java.html
		
		// initially assume all integers are prime
		
		boolean[]  isPrime;
		ArrayList<MyBigInteger>  listOfPrimes = new ArrayList<MyBigInteger>();
		
        isPrime = new boolean[n+1];
        for (int i = 2; i <= n; i++) {
            isPrime[i] = true;
        }

        // mark non-primes <= n using Sieve of Eratosthenes
        for (int factor = 2; factor*factor <= n; factor++) {

            // If factor is prime, then mark multiples of factor as nonprime
            // It suffices to consider multiples factor, factor+1, ...,  n/factor.
            if (isPrime[factor]) {
                for (int j = factor; factor*j <= n; j++) {
                    isPrime[factor*j] = false;
                }
            }
        }
		for (int i = 2;  i <= n;  i++) {
			if (isPrime[i] == true)
				listOfPrimes.add( new MyBigInteger(i,base));
		}
		return listOfPrimes;
	}
}