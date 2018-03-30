package main

import "fmt"

type matrix struct {
	data      *[][]int
	rows      int
	cols      int
	rowOffset int
	colOffset int
}

type eachElementFn func(int, int)
type transformFn func(int, int) int

/*
 * Factory method to help with creating a matrix
 */
func createMatrix(rows, cols int) matrix {
	data := make([][]int, rows)
	matrix := matrix{&data, rows, cols, 0, 0}

	for i := 0; i < rows; i++ {
		data[i] = make([]int, cols)
	}

	return matrix
}

// -- Matrix Functions -- //

/*
 * Loop over all elements of the matrix without needing
 * a loop.
 *
 * matrix.each(func(row, col int){ ... })
 *
 */
func (matrix *matrix) each(f eachElementFn) {
	for row := 0; row < matrix.rows; row++ {
		for col := 0; col < matrix.cols; col++ {
			f(row, col)
		}
	}
}

/*
 * Loop over all elements of the matrix and set
 * each element to a new value
 *
 * matrix.transform(func(row, col int) int { return ... })
 *
 */
func (matrix *matrix) transform(f transformFn) {
	for row := 0; row < matrix.rows; row++ {
		for col := 0; col < matrix.cols; col++ {
			matrix.set(row, col, f(row, col))
		}
	}
}

/*
 * Get a single row from the matrix
 *
 * row := matrix.row(0)
 *
 */
func (matrix *matrix) row(row int) *[]int {
	r := make([]int, matrix.cols)

	for col := 0; col < matrix.cols; col++ {
		r[col] = matrix.at(row, col)
	}

	return &r
}

func (matrix *matrix) set(row, col, val int) {
	data := matrix.data
	(*data)[row+matrix.rowOffset][col+matrix.colOffset] = val
}

func (matrix *matrix) at(row, col int) int {
	data := matrix.data
	return (*data)[row+matrix.rowOffset][col+matrix.colOffset]
}

// print the REAL data
func (matrix *matrix) print() {
	m := createMatrix(matrix.rows, matrix.cols)

	m.transform(func(row, col int) int {
		return matrix.at(row, col)
	})

	fmt.Println(m.data)
}

// fill the matrix with example data
func (matrix *matrix) fill() {
	matrix.transform(func(row, col int) int {
		return col + 1
	})
}

func (matrix *matrix) transpose() matrix {
	transposedMatrix := createMatrix(matrix.cols, matrix.rows)

	transposedMatrix.transform(func(row, col int) int {
		return matrix.at(col, row)
	})

	return transposedMatrix
}

/*
 * Add 2 matrices and create a new matrix
 *
 * newMatrix := matrix.add(otherMatrix)
 *
 */
func (matrix *matrix) add(b *matrix) matrix {
	c := createMatrix(matrix.rows, matrix.cols)

	c.transform(func(row, col int) int {
		return matrix.at(row, col) + b.at(row, col)
	})

	return c
}

/*
 * Multiply 2 matrices and create a new matrix
 *
 * newMatrix := matrix.multiply(otherMatrix)
 *
 */
func (matrix *matrix) multiply(b *matrix) matrix {
	c := createMatrix(matrix.rows, b.cols)
	bTrans := b.transpose()

	c.transform(func(row, col int) int {
		return dotProduct(matrix.row(row), bTrans.row(col))
	})

	return c
}

/*
 * Create a sub matrix (like a substring)
 *
 * This function is the key to keeping the memory as low as possible.
 * We give the new sub matrix a pointer to the parent matrix's data,
 * however, by setting a row and col offset along with a size,
 * we can use the other matrix methods to easily abstract the
 * fact that there is more data than needed.  Thus, sharing a pointer
 * to the data and using offsets can save tons of memory.
 *
 * subMatrix := matrix.submatrix(0, 0, 2)
 *
 */
func (matrix *matrix) submatrix(rowStart, colStart, size int) matrix {
	subMatrix := createMatrix(size, size)

	subMatrix.data = matrix.data // pointers save memory :)

	// use offsets to make referencing work like magic
	subMatrix.rowOffset = rowStart + matrix.rowOffset
	subMatrix.colOffset = colStart + matrix.colOffset

	return subMatrix
}

// -- End Matrix Functions -- //

// Split the matrix into 4 quadirants
func (matrix *matrix) quads(size int) (m11, m12, m21, m22 matrix) {
	m11 = matrix.submatrix(0, 0, size)
	m12 = matrix.submatrix(0, size, size)
	m21 = matrix.submatrix(size, 0, size)
	m22 = matrix.submatrix(size, size, size)

	return
}

// Thread wrapper for dnc on a top level quadirant (divideAndConquer)
func dncQuadThread(a, b, c, d *matrix, ch chan matrix) {
	ch <- dncQuad(a, b, c, d)
}

// divideAndConquer for a 4 quadirants
func dncQuad(a, b, c, d *matrix) matrix {
	x := dnc(a, b, false)
	y := dnc(c, d, false)

	return x.add(&y)
}

// combine 4 matrix into a single large matrix via quadirants
// [[a b]
//  [c d]]
func combineQuads(a, b, c, d *matrix, size int) matrix {
	m := createMatrix(size*2, size*2)

	a.each(func(row, col int) {
		m.set(row, col, a.at(row, col))
		m.set(row, col+size, b.at(row, col))
		m.set(row+size, col, c.at(row, col))
		m.set(row+size, col+size, d.at(row, col))
	})

	return m
}

// divide and conquer (dnc) matrix multiplication
// this is where recursion starts
func dnc(a, b *matrix, first bool) matrix {
	size := a.rows / 2

	if size == 1 {
		return a.multiply(b)
	}

	a11, a12, a21, a22 := a.quads(size)
	b11, b12, b21, b22 := b.quads(size)

	var c11, c12, c21, c22 matrix

	// if it is the first time calling dnc, start threads
	// since we are at the top most level
	if first {
		var ch1 = make(chan matrix)
		var ch2 = make(chan matrix)
		var ch3 = make(chan matrix)
		var ch4 = make(chan matrix)

		go dncQuadThread(&a11, &b11, &a12, &b21, ch1)
		go dncQuadThread(&a11, &b12, &a12, &b22, ch2)
		go dncQuadThread(&a21, &b11, &a22, &b21, ch3)
		go dncQuadThread(&a21, &b12, &a22, &b22, ch4)

		//wait
		c11 = <-ch1
		c12 = <-ch2
		c21 = <-ch3
		c22 = <-ch4

	} else {
		c11 = dncQuad(&a11, &b11, &a12, &b21)
		c12 = dncQuad(&a11, &b12, &a12, &b22)
		c21 = dncQuad(&a21, &b11, &a22, &b21)
		c22 = dncQuad(&a21, &b12, &a22, &b22)
	}

	return combineQuads(&c11, &c12, &c21, &c22, size)
}

// Dot Product of 2 Matrices
func dotProduct(a, b *[]int) int {
	c := make([]int, len(*a))
	var product int

	for i := 0; i < len(*a); i++ {
		c[i] = (*a)[i] * (*b)[i]
	}

	for i := 0; i < len(c); i++ {
		product += c[i]
	}

	return product
}

func main() {
	var size int

	// get matrix size from user
	fmt.Scanln(&size)
	fmt.Println(size)

	show := size < 16

	a := createMatrix(size, size)
	a.fill()

	b := createMatrix(size, size)
	b.fill()

	if show {
		a.print()
		b.print()
	}

	// only thread if the size is greater than 32
	c := dnc(&a, &b, size > 32)

	if show {
		c.print()
	}
}
