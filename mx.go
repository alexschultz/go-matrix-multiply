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
		return col + 1
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

// divide and conquer
func dnc(a, b *matrix) matrix {
	size := a.rows / 2
	if size == 1 {
		return a.multiply(b)
	}

	a11 := a.submatrix(0, 0, size)
	a12 := a.submatrix(0, size, size)
	a21 := a.submatrix(size, 0, size)
	a22 := a.submatrix(size, size, size)
	b11 := b.submatrix(0, 0, size)
	b12 := b.submatrix(0, size, size)
	b21 := b.submatrix(size, 0, size)
	b22 := b.submatrix(size, size, size)

	c11a := dnc(&a11, &b11)
	c11b := dnc(&a12, &b21)
	c12a := dnc(&a11, &b12)
	c12b := dnc(&a12, &b22)
	c21a := dnc(&a21, &b11)
	c21b := dnc(&a22, &b21)
	c22a := dnc(&a21, &b12)
	c22b := dnc(&a22, &b22)

	c11 := c11a.add(&c11b)
	c12 := c12a.add(&c12b)
	c21 := c21a.add(&c21b)
	c22 := c22a.add(&c22b)

	c := createMatrix(a.rows, a.cols)

	c11.each(func(row, col int) {
		c.set(row, col, c11.at(row, col))
		c.set(row, col+size, c12.at(row, col))
		c.set(row+size, col, c21.at(row, col))
		c.set(row+size, col+size, c22.at(row, col))
	})

	return c
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
	b := a
	a.fill()
	b.fill()

	dnc(&a, &b)

	// c.print()

	// var c = make(chan *[][]int)

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
