package main

import "fmt"

type matrix struct {
	data      [][]int
	rows      int
	cols      int
	rowOffset int
	colOffset int
}

type eachElementFn func(int, int)
type transformFn func(int, int) int

func createMatrix(rows, cols int) matrix {
	matrix := matrix{make([][]int, rows), rows, cols, 0, 0}

	for i := 0; i < rows; i++ {
		matrix.data[i] = make([]int, cols)
	}

	return matrix
}

func (matrix *matrix) each(f eachElementFn) {
	for row := 0; row < matrix.rows; row++ {
		for col := 0; col < matrix.cols; col++ {
			f(row, col)
		}
	}
}

func (matrix *matrix) transform(f transformFn) {
	for row := 0; row < matrix.rows; row++ {
		for col := 0; col < matrix.cols; col++ {
			matrix.set(row, col, f(row, col))
		}
	}
}

func (matrix *matrix) row(row int) *[]int {
	r := make([]int, matrix.cols)

	for col := 0; col < matrix.cols; col++ {
		r[col] = matrix.at(row, col)
	}

	return &r
}

func (matrix *matrix) set(row, col, val int) {
	matrix.data[row+matrix.rowOffset][col+matrix.colOffset] = val
}

func (matrix *matrix) at(row, col int) int {
	return matrix.data[row+matrix.rowOffset][col+matrix.colOffset]
}

func (matrix *matrix) print() {
	m := createMatrix(matrix.rows, matrix.cols)

	m.transform(func(row, col int) int {
		return matrix.at(row, col)
	})

	fmt.Println(m.data)
}

func (matrix *matrix) fill() {
	matrix.transform(func(row, col int) int {
		return 2
	})
}

func (matrix *matrix) copy() matrix {
	matrixCopy := *matrix

	return matrixCopy
}

func (matrix *matrix) transpose() matrix {
	transposedMatrix := createMatrix(matrix.cols, matrix.rows)

	transposedMatrix.transform(func(row, col int) int {
		return matrix.at(col, row)
	})

	return transposedMatrix
}

func (matrix *matrix) add(b *matrix) matrix {
	c := createMatrix(matrix.rows, matrix.cols)

	c.transform(func(row, col int) int {
		return matrix.at(row, col) + b.at(row, col)
	})

	return c
}

func (matrix *matrix) multiply(b *matrix) matrix {
	c := createMatrix(matrix.rows, b.cols)
	bTrans := b.transpose()

	c.transform(func(row, col int) int {
		return dotProduct(matrix.row(row), bTrans.row(col))
	})

	return c
}

func (matrix *matrix) submatrix(rowStart, colStart, size int) matrix {
	subMatrix := createMatrix(size, size)

	subMatrix.data = matrix.data // pointers save memory :)

	// use offsets to make referencing work like magic
	subMatrix.rowOffset = rowStart + matrix.rowOffset
	subMatrix.colOffset = colStart + matrix.colOffset

	return subMatrix
}

func (matrix *matrix) quads(size int) (m11, m12, m21, m22 matrix) {
	m11 = matrix.submatrix(0, 0, size)
	m12 = matrix.submatrix(0, size, size)
	m21 = matrix.submatrix(size, 0, size)
	m22 = matrix.submatrix(size, size, size)

	return
}

func dncQuadThread(a, b, c, d *matrix, ch chan matrix) {
	ch <- dncQuad(a, b, c, d)
}

func dncQuad(a, b, c, d *matrix) matrix {
	x := dnc(a, b, false)
	y := dnc(c, d, false)

	return x.add(&y)
}

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

// divide and conquer
func dnc(a, b *matrix, first bool) matrix {
	size := a.rows / 2

	if size == 1 {
		return a.multiply(b)
	}

	a11, a12, a21, a22 := a.quads(size)
	b11, b12, b21, b22 := b.quads(size)

	var c11, c12, c21, c22 matrix
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

	fmt.Scanln(&size)
	fmt.Println(size)

	a := createMatrix(size, size)
	a.fill()

	dnc(&a, &a, size > 32)

	// c.print()

	// var c = make(chan *[][]in)

	// a := makeMatrix(size, size)

	// for row := 0; row < size; row++ {
	// 	for col := 0; col < size; col++ {
	// 		a[row][col] = 2 * (row + col + 1)
	// 	}
	// }

	// dnc(&a, &a, c)

	// fmt.Println(<-c)

	// fmt.Println(a)
	// fmt.Println(b)
	// fmt.Println(*done)
	// fmt.Println(transpose(&b))
	// fmt.Println(*dnc(&a, &b))
	// fmt.Println(add(&a, &b))
}
