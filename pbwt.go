package main

import (
	"bufio"
	"flag"
	"os"
	"runtime/pprof"
	"time"
)

// type block struct {
// 	i int   // begin column
// 	j int   // end column
// 	k []int // row indexes
// }

type blockset struct {
	i int
	j int
	k map[int]struct{}
}

// func min(x, y int) int {
// 	if x < y {
// 		return x
// 	}
// 	return y
// }

func max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

var alphabet,
	minBlockWidth,
	minBlockRows,
	maxIndex,
	index,
	percentage,
	blocksum int

func main() {

	// flags
	alphabetPtr := flag.Int("a", 2, "Set alphabet cardinality, wildcard excluded. {0,1,*}")
	minBlockWidthPtr := flag.Int("min_w", 2, "Set minimum SNP width for blocks")
	minBlockRowsPtr := flag.Int("min_r", 2, "Set minimum rows for blocks")
	//	wildcardPtr := flag.Bool("w", true, "Set usage of wildcards in the input file")
	//	logPtr := flag.Bool("log", false, "set log boolean")
	profPtr := flag.Bool("p", false, "run with profiler")
	indexPtr := flag.Int("q", -1, "set query index to stop computation at")

	flag.Parse()

	if *profPtr {
		f, _ := os.Create("cpuprof.prof")
		defer f.Close()
		pprof.StartCPUProfile(f)
	}

	alphabet = *alphabetPtr
	minBlockWidth = *minBlockWidthPtr
	minBlockRows = *minBlockRowsPtr
	inputFile := flag.Args()[0]
	maxIndex = *indexPtr

	flag.Parse()

	// Opening file

	file, err := os.Open(inputFile)

	if err != nil {
		print("File error during LOADING (check path), aborting")
		return
	}

	info, _ := file.Stat()
	max := int(info.Size())
	// big buffer for input file/ further testing and comprehension needed
	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, max)
	scanner.Buffer(buf, max)

	scanner.Split(bufio.ScanLines)
	haplos := make([][]byte, 0)
	i := 0
	for scanner.Scan() {
		haplos = append(haplos, make([]byte, 0))
		haplos[i] = append(haplos[i], scanner.Bytes()...)
		i++
	}

	err = file.Close()
	if err != nil {
		print("File error during PARSING, aborting")
	}
	columns := len(haplos[0])
	rows := len(haplos)

	// INPUT INFO
	println(" ======================= ")
	println("| pBWT wild-blocks tool |")
	println(" ======================= ")
	println("Input", rows, "rows (samples) x", columns, "columns (SNPs)")
	println("Alphabet size:", alphabet, "\t\tMin block:", minBlockRows, "x", minBlockWidth, "\tQuery index:", index) // in further implementation the program could recognize itself input type

	// Data Init

	index = 0
	if maxIndex < 0 {
		maxIndex = columns
	}

	// // initArrays ak0 dk0 v
	ak0 := make([]int, rows)
	for i := 0; i < rows; i++ {
		ak0[i] = i
	}

	dk0 := make([]int, rows) // automaticamente tutti i valori sono inizializzati a 0

	v := make([][]int8, alphabet)
	// initArrays ak0 dk0 v

	blocksum = 0

	var since time.Duration
	var start = time.Now()

	for index < maxIndex {

		ak0, dk0, v = computeNextArrays(ak0, dk0, index, haplos)
		// ak0, dk0 = smartcollapse(ak0, dk0) //devo collassare i singoli alle li
		blockatk := computeEndingBlocks(ak0, dk0, index, v)
		blocks := len(blockatk)
		blocksum += blocks
		// if len(blockatk) > 0 {
		// 	fmt.Print("[", pivot, "]>", len(blockatk), "> ")
		// 	for i := range blockatk {
		// 		fmt.Print("<b:", blockatk[i].i, "  len:", len(blockatk[i].k), "> ;")
		// 	}
		// 	fmt.Println()
		// }

		percentage = int(float64(index) / float64(maxIndex) * 100)
		print(" ", percentage, "%  \r \t\t blocks:", blocksum, "\r")
		index++
	}

	since = time.Since(start)
	println("Started at : ", start, "\nRAN in ss: ", since)

	println("Last Arrays at index", index, ":\nak", ak0, "\ndk0", dk0, "\nv", v)

}

func computeNextArrays(ak, dk []int, k int, matrix [][]byte) ([]int, []int, [][]int8) {

	dim := len(ak)

	a := make([][]int, alphabet)
	for i := range a {
		a[i] = make([]int, 0)
	}

	d := make([][]int, alphabet)
	for i := range d {
		d[i] = make([]int, 0)
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
		a[i], d[i] = collapse(a[i], d[i])
		newdim += len(a[i])

		// a[i], d[i] = smartcollapse(a[i], d[i])

	}
	var akk []int
	for i := 0; i < alphabet; i++ {

		akk = append(akk, a[i]...)
	}
	var dkk []int
	for i := 0; i < alphabet; i++ {
		dkk = append(dkk, d[i]...)
	}

	dim = len(akk)
	v := make([][]int8, alphabet)
	for i := range v {
		v[i] = make([]int8, 1, dim)
	}

	for t := 0; t < alphabet; t++ {
		v[t][0] = 1
	}
	// se sono arrivato all'ultimo indice
	if k == len(matrix[0])-1 {
		for t := 0; t < alphabet; t++ {
			for i := 1; i < dim; i++ {
				v[t] = append(v[t], 1)
			}
		}

	} else {
		var allele int
		var prec int
		prec = int(matrix[akk[0]][k+1] - 48) //
		for i := 1; i < dim; i++ {
			allele = int(matrix[akk[i]][k+1] - 48)

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
	}

	return akk, dkk, v

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
	endingBlocks := make([]blockset, 0, len(d))

	if pivot != 0 {
		start := -1
		maxI := 0
		for i := 0; i < len(d)-1; i++ {
			if d[i] > pivot-minBlockWidth {
				if i-start >= minBlockRows {
					//check maximalitiy of maxI, pivot, start, i a
					if ends(maxI, pivot, start, i, a, v) {
						print("b2k-info: ", maxI, "  ", pivot, "  ", start, "  ", i, "\n"+
							"")
						resset := makeblockset(maxI, pivot, start, i, a)
						endingBlocks = append(endingBlocks, resset)
					}
				}
				start = i
				maxI = 0
			} else {
				maxI = max(maxI, d[i])
			}
		}
		if len(d)-1-start >= minBlockRows {
			if ends(maxI, pivot, start, len(d)-1, a, v) {

				print("blk-info: ", maxI, "  ", pivot, "  ", start, "  ", len(d)-1, "\n"+
					"")
				resset := makeblockset(maxI, pivot, start, len(d)-1, a)
				endingBlocks = append(endingBlocks, resset)
			}
		}
	}
	return endingBlocks
}

func ends(c, b, p, q int, a []int, v [][]int8) bool {
	open := false
	for t := range v {
		pivot := p + 1
		for pivot <= q {
			if v[t][pivot] != 0 {
				break
			}
			pivot++
		}
		if pivot == q+1 {
			// se scorro tutto v[t][pivot] allora il blocco e' open, posso saltarlo e andare all'index successivo
			open = true
			break
		}
	}
	return !open
}

// func computeEndingBlocks(a, d []int, pivot int, v [][]int8) []block {
// 	endingblocks := make([]block, 0, len(d))
// 	//compute openblocks
// 	for i := 1; i < len(d); i++ {
// 		x := i
// 		y := i
// 		if d[i] <= d[0]-minBlockWidth {
// 			for x > 0 && d[x] <= d[i] {
// 				x--
// 			}
// 			for y < len(d)-1 && d[y] <= d[i] {
// 				y++
// 			}
// 			y--
// 			if y-x+1 < minBlockRows {
// 				continue // se il blocco e' troppo piccolo, passo all'iterazione successiva
// 			}
// 			//check maximality with v
// 			// NUMERO DI ALLELI
// 			open := false
// 			for t := range v {
// 				// fmt.Println("DEBUG:: t", t)
// 				pivot := x + 1
// 				for pivot <= y {
// 					if v[t][pivot] != 0 {
// 						break
// 					}
// 					pivot++
// 				}
// 				if pivot == y+1 {
// 					// se scorro tutto v[t][pivot] allora il blocco e' open, posso saltarlo e andare all'index successivo
// 					open = true
// 					break
// 				}
// 			}
// 			if !open {
// 				// QUI FACCIO L'APPEND DI  block(d[i],d[0], X, Y)  ai blocchi chiusi
// 				// fmt.Println("DEBUG:: ", "d[i]=", d[i], "  d[0]=", d[0], "  x=", x, "  y=", y, "  d[x]=", d[x], "  d[y]=", d[y], "  len(d)=", len(d))
// 				resblock := makeblock(d[i], d[0], x, y, a)
// 				if len(resblock.k) >= minBlockWidth {
// 					endingblocks = append(endingblocks, resblock)
// 				}
// 			}
// 		}
// 	}
// 	return endingblocks
// }

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

// func makeblock(c, b, p, q int, a []int) block {
// 	var bl block
// 	bl.i = c
// 	bl.j = b

// 	arr := make([]int, q-p+1)
// 	for i := 0; i <= q-p; i++ {
// 		arr[i] = a[i+p]
// 	}
// 	sort.Ints(arr)
// 	arr = removeDuplicateInt(arr)
// 	bl.k = arr
// 	return bl
// }

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
