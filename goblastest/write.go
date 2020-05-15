package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"reflect"
	"strings"
)

var units map[int]*os.File
var readers map[int]*bufio.Reader
var writers map[int]*bufio.Writer

func init() {
	units = map[int]*os.File{}
	units[6] = os.Stdout
	readers = map[int]*bufio.Reader{}
	writers = map[int]*bufio.Writer{}
}

func rewind(unit int) {
	units[unit].Seek(0, 0)
}

func write(unit int, format []byte, a ...interface{}) {
again:
	for i := 1; i < len(a); i++ {
		b1, ok1 := a[i].([]byte)
		for j := i + 1; j < len(a); j++ {
			b2, ok2 := a[j].([]byte)
			if ok1 && ok2 {
				b1 = append(b1, b2...)
				a = append(a[:i], append([]interface{}{b1}, a[i+1:]...)...)
				goto again
			}
		}
	}

	for i := range a {
		if str, ok := a[i].([]byte); ok {
			a[i] = string(str)
		}
	}

	fmt.Fprintf(units[unit], string(format), a...)
}

func writeString(unit int, format string, a ...interface{}) {
again:
	for i := 1; i < len(a); i++ {
		b1, ok1 := a[i].([]byte)
		for j := i + 1; j < len(a); j++ {
			b2, ok2 := a[j].([]byte)
			if ok1 && ok2 {
				b1 = append(b1, b2...)
				a = append(a[:i], append([]interface{}{b1}, a[i+1:]...)...)
				goto again
			}
		}
	}

	for i := range a {
		if str, ok := a[i].([]byte); ok {
			a[i] = string(str)
		}
	}

	fmt.Fprintf(units[unit], format, a...)
}

func openWrite(unit int, file string) {
	f, err := os.Create(file)
	if err != nil {
		panic(err)
	}

	units[unit] = f
}

func openRead(unit int, file string) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}

	units[unit] = f
	readers[unit] = bufio.NewReader(f)
}

func _close(unit int) {
	units[unit].Close()
	delete(units, unit)
}

func read(unit int, format string, a ...interface{}) {

	format = strings.TrimSpace(strings.ToLower(format))

	// Change from %15.4f to %15f
	for i := 1; i < 20; i++ {
		for j := 1; j <= i; j++ {
			format = strings.Replace(format,
				fmt.Sprintf("%c%d.%df", '%', i, j),
				fmt.Sprintf("%cf", '%'),
				-1)
		}
	}

	_, err := fmt.Fscanf(units[unit], format, a)
	if err != nil {
		var types string
		for i := range a {
			types += fmt.Sprintf("|%s|", reflect.TypeOf(a[i]))
		}
		panic(fmt.Errorf("READ error for format `%s` : %v\nValues = %v",
			format, err, types))
	}
}

func readln(r *bufio.Reader) string {
	line, _, err := r.ReadLine()
	if err != nil && err != io.EOF {
		panic("")
	} else if err == io.EOF {
		return ""
	}
	return string(line)
}

func readval(r *bufio.Reader) string {
	line, err := r.ReadString(' ')
	if err != nil {
		panic("")
	}
	return line
}
