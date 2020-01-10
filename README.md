# Big-Integer

An implementation of a Big Integer class, whcih includes several implementations of grade school arithmetic, base conversion, and prime factorization algorithms on big integers. The three main methods are:

dividedBy

    An implementation of long division. The time complexity is O(N^2).
    
convert

    Converts a Big Integer from one base to another. The convert method is called by a MyBigInteger object that
    represents a positive integer in some base, and returns the same positive integer represented in the new base.
    The bases can be any of {2, 3, 4, …, 10}. 

primeFactors

    The method returns an ArrayList<MyBigInteger> whose elements are the prime factors of
    𝑚, where 𝑚 is ‘this’ MyBigInteger object. The number of copies of the prime factor in the
    returned list is the order of that prime factor, and the prime factors in the list goes from
    smallest to largest. For example, if 𝑚 = 24, the method returns a list (2, 2, 2, 3). If 𝑚 is
    prime or if the method cannot find any prime factors, then the returned list will just contain one
    element, namely, the returned list will be (𝑚).
