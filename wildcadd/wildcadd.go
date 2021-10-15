package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"time"
)

var Reader io.Reader
var wildrate float64
var outfolder string

func main() {

	rate := flag.Float64("r", 0, "set wildcard insertion rate, double")
	folder := flag.String("o", "./", "choose the output folder")
	flag.Parse()

	wildrate = *rate
	outfolder = *folder

	inputFile := flag.Args()[0]
	file, err := os.Open(inputFile)
	if err != nil {
		// add error message
		return
	}

	info, _ := file.Stat()
	max := int(info.Size())

	// count columns
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
		log.Fatalf("File error, aborting")
	}
	rows := i
	columns := len(haplos[0])

	fmt.Println("INPUT------")
	fmt.Println("n of rows:", rows)
	fmt.Println("n of col:", columns)

	// end count columns
	count := (rows * columns * int(wildrate*100)) / 100

	// gen the rand position in count
	positions := make(map[int]struct{})
	fmt.Println("n wild:", count)
	for x := 0; x < count; x++ {
		rand.Seed(time.Now().UnixNano())
		positions[rand.Intn(rows*columns)] = struct{}{}
	}
	fmt.Println("-----Inserting Wildcards")
	for x := range positions {
		j := x / rows
		i := x % rows
		haplos[i][j] = 42
	}

	name := outfolder + "wild_" + fmt.Sprint(wildrate) + "_" + info.Name()
	f, _ := os.Create(name)
	fmt.Println("-----Writing in file: ", f.Name())
	defer f.Close()
	w := bufio.NewWriterSize(f, max)

	for i := 0; i < rows; i++ {
		w.Write(haplos[i])
		if i != rows-1 {
			w.WriteString("\n")
		}
	}

	w.Flush()
	f.Sync()
}
