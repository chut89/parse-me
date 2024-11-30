package com.example.metalchemist;

import java.util.stream.IntStream;
import java.nio.charset.StandardCharsets;
import java.io.BufferedReader;
import java.io.Reader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.IOException;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

import java.util.Map;
import java.util.TreeMap;

public class LoadAndRunTest {
    private String filename;
    private String blankHeader;

    public LoadAndRunTest(String filename) {
        this.filename = filename;
    }

    // To debug a jar file run: java -Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=y,address=1044 -jar target/parse-me-1.0-SNAPSHOT.jar /home/chut/workspace/codewars/metalChemist/parse-me/testdata/testcase_20241104
    // To connect to a remote debugger on port 1044 run: jdb -attach localhost:1044
    public void load() {
        try (final BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(filename), StandardCharsets.UTF_8.name()))) {
            String line = null;
            Map<String, Integer> actual = new TreeMap<>();
            boolean passed = true;
            while ((line = in.readLine()) != null) {
                if (!line.startsWith("-") && !line.startsWith("Molecule has been set to:")) {
                    if (!line.startsWith("Should lead to: ")) {
                        actual = new ParseHer(line).parse();
                    } else {
                        System.out.println(line);
                        Matcher resultMatcher = Pattern.compile("([a-zA-Z]+):(?: )?(\\d+)").matcher(line.substring("Should lead to: ".length()));
                        while (resultMatcher.find()) {
                            if (actual.get(resultMatcher.group(1)) != Integer.parseInt(resultMatcher.group(2))) {
                                passed = false;
                                break;
                            }
                        }
                        if (passed) {
                            System.out.println("----------------------------");
                        } else {
                            System.err.println("But actual result is: " + actual);
                            throw new RuntimeException("Test failed: ");
                        }
                    }
                }
            }
    
        } catch (IOException exception) {
            System.out.println("Exception occurred on reading file name:" + filename);
        }
    }
    
    public static void main(String[] args) {
        new LoadAndRunTest(args[0]).load();
    }
}
