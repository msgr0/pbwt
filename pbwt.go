package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime/pprof"

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

		if (allele <= 9) && (allele >= 0) {

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
		// a[i], d[i] = collapse(a[i], d[i])
		// smartcollapse(a[i], d[i])

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
func smartcollapse(ak, dk []int) ([]int, []int) {
	a := make([]int, len(ak))
	d := make([]int, len(dk))
	copy(a, ak)
	copy(d, dk)

	if len(a) > 1 {
		i := 0
		j := 0
		for i < len(a)-1 { // scorro tutto l'array a
			// se trovo due valori diversi passo al successivo
			if a[i] != a[i+1] && i <= j+1 {
				i++
				j = i
			} else if a[i] != a[i+1] && i > j+1 {
				a = append(a[:j+1], a[i:]...)
				d = append(d[:j+1], d[i:]...)
				i++
				j = i
			} else { // mi trovo davanti alla prima occorrenza uguale append a[:j], a[i:] if i > j
				i++
			}
		}
	}
	return a, d
}
func collapse(a, d []int) ([]int, []int) {
	if len(a) > 1 {
		ac := []int{}
		dc := []int{}
		// ac := make([]int, 0, len(a))
		// dc := make([]int, 0, len(d))

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
	return a, d
}

func computeEndingBlocks(a, d []int, pivot int, v [][]int8) []blockset {
	endingblocks := make([]blockset, 0, len(d))
	//compute openblocks
	for i := 1; i < len(d); i++ {
		x := i
		y := i

		if d[i] <= d[0]-minBlockWidth {
			for x > 0 && d[x] <= d[i] {
				x--
			}

			for y < len(d)-1 && d[y] <= d[i] {
				y++
			}
			y--

			if y-x+1 < minBlockRows {
				continue // se il blocco e' troppo piccolo, passo all'iterazione successiva
			}

			//check maximality with v
			// NUMERO DI ALLELI
			open := false
			for t := range v {
				// fmt.Println("DEBUG:: t", t)
				pivot := x + 1
				for pivot <= y {
					if v[t][pivot] != 0 {
						break
					}
					pivot++
				}
				if pivot == y+1 {
					// se scorro tutto v[t][pivot] allora il blocco e' open, posso saltarlo e andare all'index successivo
					open = true
					break
				}

			}
			if !open {
				// QUI FACCIO L'APPEND DI  block(d[i],d[0], X, Y)  ai blocchi chiusi
				// fmt.Println("DEBUG:: ", "d[i]=", d[i], "  d[0]=", d[0], "  x=", x, "  y=", y, "  d[x]=", d[x], "  d[y]=", d[y], "  len(d)=", len(d))
				// resblock := makeblock(d[i], d[0], x, y, a)
				// if len(resblock.k) >= minBlockWidth {
				// 	endingblocks = append(endingblocks, resblock)
				// }
				resset := makeblockset(d[i], d[0], x, y, a)
				if len(resset.k) >= minBlockWidth {
					endingblocks = append(endingblocks, resset)
				}

			}
		}

	}

	return endingblocks
}

func makeblock(c, b, p, q int, a []int) block {
	var bl block
	bl.i = c
	bl.j = b

	arr := make([]int, q-p+1)
	for i := 0; i <= q-p; i++ {
		arr[i] = a[i+p]
	}
	// sort.Ints(arr)
	// arr = removeDuplicateInt(arr)
	bl.k = arr

	return bl

}

func makeblockset(c, b, p, q int, a []int) blockset {
	var bls blockset
	bls.i = c
	bls.j = b

	set := make(map[int]struct{})

	for i := 0; i <= q-p; i++ {
		set[a[i+p]] = struct{}{}
	}

	bls.k = set
	return bls

}

// func removeDuplicateInt(intSlice []int) []int {
// 	allKeys := make(map[int]bool)
// 	list := []int{}
// 	for _, item := range intSlice {
// 		if _, value := allKeys[item]; !value {
// 			allKeys[item] = true
// 			list = append(list, item)
// 		}
// 	}
// 	return list
// }

type block struct {
	i int   // begin column
	j int   // end column
	k []int // row indexes
}

type blockset struct {
	i int
	j int
	k map[int]struct{}
}

var alphabet int
var wildcard bool
var minBlockWidth int
var minBlockRows int
var index int

func main() {

	// flags
	alphabetPtr := flag.Int("a", 2, "Set alphabet cardinality, wildcard excluded. {0,1,*}")
	minBlockWidthPtr := flag.Int("min_w", 2, "Set minimum SNP width for blocks")
	minBlockRowsPtr := flag.Int("min_r", 2, "Set minimum rows for blocks")
	wildcardPtr := flag.Bool("w", true, "Set usage of wildcards in the input file")
	logPtr := flag.Bool("log", false, "set log boolean")
	profPtr := flag.Bool("p", false, "run with profiler")
	indexPtr := flag.Int("q", -1, "set query index to stop computation at")

	flag.Parse()

	if !*logPtr {
		log.SetLevel(0)
	}

	if *profPtr {
		f, _ := os.Create("cpuprof.prof")
		defer f.Close()
		pprof.StartCPUProfile(f)
	}

	alphabet = *alphabetPtr
	wildcard = *wildcardPtr
	minBlockWidth = *minBlockWidthPtr
	minBlockRows = *minBlockRowsPtr
	inputFile := flag.Args()[0]
	index = *indexPtr

	// serve questo flag parse?

	flag.Parse()

	// Opening file

	file, err := os.Open(inputFile)

	if err != nil {
		return
	}
	info, _ := file.Stat()
	max := int(info.Size())
	// big buffer for input file/ further testing and comprehension needed
	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, max)
	scanner.Buffer(buf, max)

	scanner.Split(bufio.ScanLines)
	var haplos []string
	for scanner.Scan() {
		haplos = append(haplos, scanner.Text())
	}

	err = file.Close()
	if err != nil {
		log.Fatalf("File error, aborting")
	}
	columns := len(haplos[0][:])
	rows := len(haplos)

	// INPUT INFO

	fmt.Println("========================================================")
	fmt.Println("||                 pBWT wild-blocks tool              ||")
	fmt.Println("========================================================")
	fmt.Println("Input matrix is", rows, "rows (samples) x", columns, "columns (SNPs)")
	fmt.Println("Assuming: alphabet size ==", alphabet) // in further implementation the program could recognize itself input type
	fmt.Println("Assuming: wildcard presence ==", wildcard)
	fmt.Println("Assuming: min block size ==", minBlockRows, "rows x", minBlockWidth, "columns")
	fmt.Println("Assuming: queryIndex ==", index)

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
	var maxarray int
	if index < 0 {
		maxarray = len(haplos[0][:])
	} else {
		maxarray = index
	}

	var since time.Duration
	var start = time.Now()
	blocksum := 0
	for pivot < maxarray {

		ak0, dk0 = computeNextArrays(ak0, dk0, pivot, haplos)
		ak0, dk0 = smartcollapse(ak0, dk0) //devo collassare i singoli alleli
		// CALCOLO DEL BITVECTOR
		v = computeBitVectors(ak0, dk0, pivot, haplos)
		// CALCOLO DEI BLOCCHI TERMINANTI IN PIVOT
		blockatk := computeEndingBlocks(ak0, dk0, pivot, v)

		pivot++
		blocks := len(blockatk)
		blocksum += blocks
		// if blocks > 0 {
		// 	fmt.Print("[", pivot, "]>", blocks, " ")
		// 	for i := range blockatk {
		// 		fmt.Print("<b:", blockatk[i].i, "  len:", len(blockatk[i].k), "> ;")
		// 	}
		// 	fmt.Println()
		// }
	}

	since = time.Since(start)
	fmt.Println("ENDED with a total of ", blocksum)
	fmt.Println("Started at : ", start, "\nRAN in ss: ", since)
	// fmt.Println("Last Arrays\nak", ak0, "\ndk0", dk0, "\nv", v)
	if *profPtr {
		pprof.StopCPUProfile()
	}

}
