package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"

	// "sort"
	"time"

	log "github.com/sirupsen/logrus"
)

func computeBitVectors(ak, dk []int, k int, matrix []string) [][]int8 {
	dim := len(ak)

	v := make([][]int8, alphabet)
	for i := range v {
		v[i] = make([]int8, 1, dim)
	}
	// init first position
	for t := 0; t < alphabet; t++ {
		v[t][0] = 1
	}

	if k == len(matrix[0])-1 {
		for t := 0; t < alphabet; t++ {
			for i := 1; i < dim; i++ {
				v[t] = append(v[t], 1)
			}
		}
		return v
	}

	var allele int
	var prec int
	prec = int(matrix[ak[0]][k+1] - 48) //
	log.Println("Prec ", prec)
	for i := 1; i < dim; i++ {
		allele = int(matrix[ak[i]][k+1] - 48)
		// if allele == prec {
		// 	for t := range v {
		// 		v[t] = append(v[t], 0)
		// 	}

		// }

		for t := 0; t < alphabet; t++ {
			if allele == prec {
				v[t] = append(v[t], 0)
			} else if allele > 9 {
				if t == prec {
					v[t] = append(v[t], 0)
				} else {
					v[t] = append(v[t], 1)
				}
			} else {
				v[t] = append(v[t], 1)
			}
		}
		prec = allele
	}
	return v
}

func computeNextArrays(ak, dk []int, k int, matrix []string) ([]int, []int) {
	dim := len(ak)
	// allocing dim size, that's not really memory wise to be honest
	// go slices got me covered ehehehe
	a := make([][]int, alphabet)
	for i := range a {
		a[i] = make([]int, 0, dim)
	}

	d := make([][]int, alphabet)
	for i := range d {
		d[i] = make([]int, 0, dim)
	}

	p := make([]int, alphabet)
	u := make([]int, alphabet)

	for i := 0; i < alphabet; i++ {
		u[i] = 0
		p[i] = k + 1
	}

	var allele int
	for i := 0; i < dim; i++ {
		allele = int(matrix[ak[i]][k] - 48)

		for l := 0; l < alphabet; l++ {
			if dk[i] > p[l] {
				p[l] = dk[i]
			}
		}

		if (allele < 10) && (allele >= 0) {

			a[allele] = append(a[allele], ak[i])
			d[allele] = append(d[allele], p[allele])
			p[allele] = 0
			u[allele] = u[allele] + 1

		} else {
			for m := 0; m < alphabet; m++ {
				a[m] = append(a[m], ak[i])
				d[m] = append(d[m], p[m])
				p[m] = 0
				u[m] = u[m] + 1
			}
		}

	}

	newdim := 0
	for i := 0; i < alphabet; i++ {
		newdim += len(a[i])
	}
	var akk []int
	for i := 0; i < alphabet; i++ {
		akk = append(akk, a[i]...)
	}
	var dkk []int
	for i := 0; i < alphabet; i++ {
		dkk = append(dkk, d[i]...)
	}

	return akk, dkk

}

func collapse(a, d []int) ([]int, []int) {
	ac := make([]int, 0, len(a))
	dc := make([]int, 0, len(d))

	j := 0
	pivot := 0
	for j < len(a)-1 {
		if a[j] == a[j+1] {
			j++
		} else {
			ac = append(ac, a[pivot])
			dc = append(dc, d[pivot])
			j++
			pivot = j
		}

	}

	ac = append(ac, a[pivot])
	dc = append(dc, d[pivot])

	return ac, dc
}

func computeBlocks(a, d []int, v [][]int8) []block {
	// scorro gli indici,
	// cerco la finestra per un blocco
	// controllo per quegli indici se bit vector chiude il blocco
	// se chiude appendo ai blocchi chiusi in output
	// se non chiude cambio indice

	return
}

func availableBlocks(a, d []int) []block {
	availblock := make([]block, 0, len(d))
	for i := 1; i < len(d)-1; i++ {
		x := i - 1
		y := i + 1
		if d[i] <= d[0]-minBlockWidth {
			// continuo ad espandere i blocchi con row sopra e sotto in base al livello di D e ai bound
			for x > 0 && d[x] <= d[i] {
				x--
			}

			for y < len(d) && d[y] <= d[i] {
				y++
			}

			// controllo se si puo espandere
			// se si lascio perdere
			// se no chiudo il blocco e lo appendo in uscita

			// block trim
			newblock := makeblock(d[i], d[0], a[x:y])
			log.Println("appending block: ", newblock)
			availblock = append(availblock, newblock)

		}

	}

	return availblock

	//TODO
	//TRIM to 1 instance of block
	// check maximality with BItVector

}

func computeEndingBlocks(a, d []int, pivot int, v [][]int8) []block {
	endingblocks := make([]block, 0, len(d))

	return endingblocks
}

type block struct {
	i int   // begin column
	j int   // end column
	k []int // row indexes
}

type openblock struct {
	a int // column index beg
	b int // column index end
	i int // row index begin
	j int // row index end
}

var alphabet int
var wildcard bool
var minBlockWidth int
var minBlockRows int

func main() {
	// flags
	alphabetPtr := flag.Int("a", 2, "Set alphabet cardinality, wildcard excluded. {0,1,*}")
	minBlockWidthPtr := flag.Int("min_w", 2, "Set minimum SNP width for blocks")
	minBlockRowsPtr := flag.Int("min_r", 2, "Set minimum rows for blocks")
	wildcardPtr := flag.Bool("w", true, "Set usage of wildcards in the input file")
	logPtr := flag.Bool("log", false, "set log boolean")

	flag.Parse()

	if !*logPtr {
		log.SetLevel(0)
	}

	alphabet = *alphabetPtr
	wildcard = *wildcardPtr
	minBlockWidth = *minBlockWidthPtr
	minBlockRows = *minBlockRowsPtr
	inputFile := flag.Args()[0]

	// serve questo flag parse?

	flag.Parse()

	// Opening file

	file, err := os.Open(inputFile)

	if err != nil {
		return
	}

	// big buffer for input file/ further testing and comprehension needed
	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, 64*5000)
	scanner.Buffer(buf, 1024*1024)

	scanner.Split(bufio.ScanLines)
	var haplos []string
	for scanner.Scan() {
		haplos = append(haplos, scanner.Text())
	}

	err = file.Close()
	if err != nil {
		log.Fatalf("File error, aborting")
	}
	columns := len(haplos[:])
	rows := len(haplos)

	// INPUT INFO

	fmt.Println("========================================================")
	fmt.Println("||                 pBWT wild-blocks tool              ||")
	fmt.Println("========================================================")
	fmt.Println("Input matrix is", rows, "rows (samples) x", columns, "columns (SNPs)")
	fmt.Println("Assuming: alphabet size ==", alphabet) // in further implementation the program could recognize itself input type
	fmt.Println("Assuming: wildcard presence ==", wildcard)
	fmt.Println("Assuming: min block size ==", minBlockRows, "rows x", minBlockWidth, "columns")

	// initArrays ak0 dk0 v

	ak0 := make([]int, 0, rows)
	for i := 0; i < rows; i++ {
		ak0 = append(ak0, i)
	}

	dk0 := make([]int, 0, rows)
	for i := 0; i < rows; i++ {
		dk0 = append(dk0, 0)
	}

	v := make([][]int8, alphabet)

	pivot := 0
	maxarray := len(haplos[:])

	var since time.Duration
	var start = time.Now()

	for pivot < maxarray {

		ak0, dk0 = computeNextArrays(ak0, dk0, pivot, haplos)
		ak0, dk0 = collapse(ak0, dk0)
		// CALCOLO DEL BITVECTOR
		v = computeBitVectors(ak0, dk0, pivot, haplos)
		// CALCOLO DEI BLOCCHI TERMINANTI IN PIVOT
		blockatk := computeEndingBlocks(ak0, dk0, pivot, v)

		fmt.Println("Blocks ending at k = ", pivot, ":")

		pivot++
	}

	since = time.Since(start)
	fmt.Println("Started at : ", start, "\nRAN in ss: ", since)
	fmt.Println("Last Arrays\nak", ak0, "\ndk0", dk0, "\nv", v)

}
