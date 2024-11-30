package com.example.metalchemist;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.Optional;
import java.util.stream.IntStream;

public class ParseHer {
    
    //                            Number      :   1       2      3      4
    final static private String[] RADICALS    = {"meth", "eth", "prop", "but",   "pent",  "hex",  "hept",  "oct",  "non",  "dec",  "undec",  "dodec",  "tridec",  "tetradec",  "pentadec",  "hexadec",  "heptadec",  "octadec",  "nonadec"},
                                  MULTIPLIERS = {        "di",  "tri",  "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca", "undeca", "dodeca", "trideca", "tetradeca", "pentadeca", "hexadeca", "heptadeca", "octadeca", "nonadeca"},
                                  
                                  SUFFIXES    = {         "ol",      "al", "one", "oic acid", "carboxylic acid",                "oate",             "ether", "amide", "amine", "imine", "benzene", "thiol",    "phosphine", "arsine"},
                                  PREFIXES    = {"cyclo", "hydroxy",       "oxo",             "carboxy",         "oxycarbonyl", "oyloxy", "formyl", "oxy",   "amido", "amino", "imino", "phenyl",  "mercapto", "phosphino", "arsino", "fluoro", "chloro", "bromo", "iodo"}; // 19 in total
    
    // Note that alkanes, alkenes alkynes, and akyles aren't present in these lists
    // constant member as elements
    final private String CARBON = "C";
    final private String OXYGEN = "O";
    final private String HYDROGEN = "H";
    final private String NITROGEN = "N";

    private String molec = "";

    // member variables that represent the state of molecule computation
    private int numberOfCarbons = 0;
    private int numberOfCyclos = 0;
    private int numberOfAlkeneBounds = 0;
    private int numberOfAlkyneBounds = 0;
    private int numberOfFunctionBounds = 0;
    private Map<String, Integer> result = null;
    //private String suffix;

    // The following member variables are for esterSuffix
    private int esterSuffixNumberOfCarbons = 0;
    private int esterSuffixNumberOfCyclos = 0;
    private int esterSuffixNumberOfAlkeneBounds = 0;
    private int esterSuffixNumberOfAlkyneBounds = 0;
    private int esterSuffixNumberOfFunctionBounds = 0;
    private Map<String, Integer> esterSuffixResult = null;
    private int esterSuffixMultiplier = -1;
    
    // precompiled usefull pattern
    private Pattern alkylPattern = null;
    private Pattern functionPattern = null;

    public ParseHer(String name) {
        // Do whatever you want here...
        this.molec = name;

        System.out.println("Molecule has been set to: " + molec);

        result = new TreeMap<>();
        esterSuffixResult = new TreeMap<>();

        StringBuilder alkylPatternBuilder = new StringBuilder("(((?:\\d+,?)+-)?(?:");
        // (1-3-4-)tripent(en|yn|enyn)(oxycarbon)yl
        Arrays.stream(MULTIPLIERS).forEach(mult -> alkylPatternBuilder.append(mult + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")?(?:");
        alkylPatternBuilder.append(PREFIXES[0]);
        alkylPatternBuilder.append(")?(?:");
        Arrays.stream(RADICALS).forEach(radical -> alkylPatternBuilder.append(radical + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        // (1-3-)(di)hex-2-en-3,4-triynyl
        alkylPatternBuilder.append(")(?:en|yn|enyn)?(oxycarbon)?yl|(((?:\\d+,?)+-)?(");
        Arrays.stream(MULTIPLIERS).forEach(mult -> alkylPatternBuilder.append(mult + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")?(?:");
        alkylPatternBuilder.append(PREFIXES[0]);
        alkylPatternBuilder.append(")?(?:");
        Arrays.stream(RADICALS).forEach(radical -> alkylPatternBuilder.append(radical + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")((?:-(?:\\d+,?)+-(?:");
        Arrays.stream(MULTIPLIERS).forEach(mult -> alkylPatternBuilder.append(mult + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")?(?:en|yn|enyn))+)(oxycarbon)?yl)|");
        // (1)-dodecen-3,6,8-triynyl
        alkylPatternBuilder.append("((?:\\d+,?)+-)?(?:");
        Arrays.stream(MULTIPLIERS).forEach(mult -> alkylPatternBuilder.append(mult + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")?(?:");
        alkylPatternBuilder.append(PREFIXES[0]);
        alkylPatternBuilder.append(")?(?:");
        Arrays.stream(RADICALS).forEach(radical -> alkylPatternBuilder.append(radical + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")en-(?:\\d+,?)+-(?:");
        Arrays.stream(MULTIPLIERS).forEach(mult -> alkylPatternBuilder.append(mult + "|"));
        alkylPatternBuilder.deleteCharAt(alkylPatternBuilder.length() - 1);
        alkylPatternBuilder.append(")?ynyl");
        alkylPatternBuilder.append(")-?");
        alkylPattern = Pattern.compile(alkylPatternBuilder.toString());

        // propoxycarbonyl   ethenoxycarbonyl

        StringBuilder functionPatternBuilder = new StringBuilder("(((?:(?:\\d+,?)+-)?)(?:");
        Arrays.stream(MULTIPLIERS).forEach(mult -> functionPatternBuilder.append(mult + "|"));
        functionPatternBuilder.deleteCharAt(functionPatternBuilder.length() - 1);
        functionPatternBuilder.append(")?(?:");
        Arrays.stream(PREFIXES).filter(prefix -> !prefix.equals(PREFIXES[0]) && !prefix.equals(PREFIXES[4]) && !prefix.equals(PREFIXES[5]) && !prefix.equals(PREFIXES[7]))
            .forEach(prefix -> functionPatternBuilder.append(prefix + "|")); // We have special validations for -oxycarbonyl -oyloxy and -oxy
        functionPatternBuilder.deleteCharAt(functionPatternBuilder.length() - 1);
        functionPatternBuilder.append("))-?");
        functionPattern = Pattern.compile(functionPatternBuilder.toString());
    }

    // Putting [.] inside Prefix pattern is dangerous and prone to confusion: []...[] for example
    private final String GROUP1_EXP = "((?:(?:(?:\\d+,?)+-)?[a-zA-Z]*(?:(?:(?:(?:-(?:\\d+,?)+-[a-zA-Z]*)?(?:en|yn|enyn))*|an)(?:oyloxy|oxy)|phenyl|formyl|(?:-(?:\\d+,?)+-[a-zA-Z]*(?:en|yn|enyn)|oxycarbon)*yl|bromo|imino|amino|phosphino|arsino|amido|iodo|mercapto|fluoro|chloro|oxo|carboxy|hydroxy|oxy)(?:-| )?)*)";
    // BIG TODO: rework GROUP2_EXP; suffix can also contain prefix functional group
    // Possible approach: work out bracket expression first!!!
    private final String GROUP2_EXP = "([a-zA-Z]*(?:(?:-(?:\\d+,?)+-)*[a-zA-Z]*(?:en|yn|enyn|an))*(?:(?:-(?:\\d+,?)+-)*[a-zA-Z]*(?:ene|yne)|ane|(?:-(?:\\d+,?)+-)?[a-zA-Z]*(?:ol|thiol|imine|amine|phosphine|arsine|al|one|amide|oic acid|carboxylic acid|ether|oate|benzene)))";

    public Map<String,Integer> parse() {
        // Parse the name given as argument in the constructor and output the Map representing the raw formula
        if (Arrays.stream(molec.split("[a-zA-Z0-9-\\[\\] ,]+")).count() > 0) {
            throw new RuntimeException("Token only accepts alphanumeric and [] and dash and space characters, molec: " + molec);
        }

        int startIdx = 0;
        boolean cyclic = false;

        Matcher esterSuffixMatcher = Pattern.compile("((?:\\d+,?)+-)+([a-zA-Z]*)oate").matcher(molec);
        if (esterSuffixMatcher.find()) {
            String digitToken = esterSuffixMatcher.group(1);
            int digitCounter = (int) Arrays.stream(digitToken.substring(0, digitToken.indexOf("-")).split(",")).count();
            if (digitCounter > 1 && esterSuffixMatcher.group(2).equals(MULTIPLIERS[digitCounter - 2])) {
                esterSuffixMultiplier = digitCounter;
            } else if (digitCounter > 1 && !esterSuffixMatcher.group(2).startsWith(MULTIPLIERS[digitCounter - 2])) {
                throw new RuntimeException("Something goes wrong with multiplier of 'oate' function, molec: " + molec);
            }
        }

        RelativeIndexHolder spaceIdxHolder = new RelativeIndexHolder(molec.indexOf(" "));
        StringBuilder ramificationBuilder = new StringBuilder(molec);

        while (ramificationBuilder.indexOf("[") >= 0 && ramificationBuilder.indexOf("]") >= 0) {
            int openBrCounter = 1;
            int openBrIdx = ramificationBuilder.indexOf("[");
            int charIdx = openBrIdx + 1;

            int lastDashIdx = ramificationBuilder.substring(0, openBrIdx).lastIndexOf("-");
            int spaceIdx = ramificationBuilder.indexOf(" ");
            Optional<Integer> optMultiplierIdx = null;
            int multiple = -1;
            startIdx = lastDashIdx + 1;
            if (lastDashIdx >= 0) {
                int firstDashIdx = ramificationBuilder.substring(0, lastDashIdx).lastIndexOf("-");
                int posCount = -1;
                if (firstDashIdx >= 0) {
                    posCount = (int) Arrays.stream(ramificationBuilder.substring(firstDashIdx + 1, lastDashIdx).split(",")).count();
                } else if (!Character.isDigit(ramificationBuilder.charAt(0)) && !(spaceIdx >= 0 && spaceIdx < lastDashIdx)) {
                    throw new RuntimeException(String.format("Dangle dash followed by subchain, molec: %s", molec));
                } else {
                    posCount = (int) Arrays.stream(ramificationBuilder.substring(Character.isDigit(ramificationBuilder.charAt(0)) ? 0 : spaceIdx + 1, lastDashIdx).split(",")).count();
                }
                if (ramificationBuilder.substring(lastDashIdx + 1, openBrIdx).equals(posCount > 1 ? MULTIPLIERS[posCount - 2] : "")) {
                    multiple = posCount;
                } else if (posCount > 1 && !ramificationBuilder.substring(lastDashIdx + 1).startsWith(MULTIPLIERS[posCount - 2])) {
                    throw new RuntimeException("Multipliers of subchain doesn't match the number of proceeding positions, molec:" + molec);
                } else if (posCount > 1) {
                    startIdx += MULTIPLIERS[posCount - 2].length();
                }
            }
            if (multiple < 0) {
                final String multiplierCand = ramificationBuilder.substring(startIdx, openBrIdx);
                optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                    .boxed()
                    .sorted((i1, i2) -> i2 - i1)
                    .filter(i -> multiplierCand.endsWith(MULTIPLIERS[i]))
                    .findFirst();
            
                multiple = optMultiplierIdx.orElse(-1) + 2;
            }

            while (charIdx < ramificationBuilder.length() && openBrCounter > 0) {
                if (ramificationBuilder.charAt(charIdx) == '[') {
                    ++openBrCounter;
                } else if (ramificationBuilder.charAt(charIdx) == ']') {
                    --openBrCounter;
                }
                charIdx++;
            }
            int closeBrIdx = charIdx - 1;
            RelativeIndexHolder newSpaceIdxHolder = new RelativeIndexHolder(spaceIdxHolder.getRelIdx());
            if (spaceIdxHolder.getRelIdx() >= 0) {
                newSpaceIdxHolder.updateInsideBrackets(openBrIdx, closeBrIdx);
                spaceIdxHolder.updateOnShortenedString(openBrIdx, closeBrIdx);
            }
            parseSubchain(ramificationBuilder.substring(openBrIdx + 1, closeBrIdx), multiple, newSpaceIdxHolder);
            ramificationBuilder.delete(openBrIdx, closeBrIdx + 1);
        }

        Matcher matcher = Pattern.compile(GROUP1_EXP + GROUP2_EXP).matcher(ramificationBuilder.toString());

        if (matcher.find()) {
            String matchedWholeWord = matcher.group(1);
            if (matchedWholeWord != null && !matchedWholeWord.isEmpty()) {
                // include this map into this.result
                parseSubchain(matchedWholeWord, 1, spaceIdxHolder);

                if (esterSuffixMultiplier > 1 && matchedWholeWord.contains("yl ")) {
                    numberOfCarbons += esterSuffixNumberOfCarbons * esterSuffixMultiplier;
                    numberOfCyclos += esterSuffixNumberOfCyclos * esterSuffixMultiplier;
                    numberOfAlkeneBounds += esterSuffixNumberOfAlkeneBounds * esterSuffixMultiplier;
                    numberOfAlkyneBounds += esterSuffixNumberOfAlkyneBounds * esterSuffixMultiplier;
                    numberOfFunctionBounds += esterSuffixNumberOfFunctionBounds * esterSuffixMultiplier;
                    for (Map.Entry<String, Integer> entry : esterSuffixResult.entrySet()) {
                        String origKey = entry.getKey();
                        result.put(origKey, (result.containsKey(origKey) ? result.get(origKey) : 0) + entry.getValue() * esterSuffixMultiplier);
                    }
                }
            }

            int posCount = 0;
            Optional<Integer> optRadicalIdx;
            String multiplier = null;

            Matcher suffixMatcher = Pattern.compile(
                "(?:benzene|" +
                "(([a-zA-Z]+)(en|yn|enyn|an)|([a-zA-Z]+-)(?:(?:(?:\\d+,?)+-)?(?:[a-zA-Z]{2,}?)?(?:en|yn|enyn|an)-?)*)" +
                "((?:(?:\\d+,?)+-)?([a-zA-Z]{2,}?)?(?:ene|yne|thiol|ol|imine|al|amine|phosphine|arsine|one|amide|oic acid|carboxylic acid|ether|oate)|ane)" +
                ")"
            ).matcher(matcher.group(2));

            String pretoken = null;
            String suftoken = null;

            String radicalInMainChainFull = null;

            if (!suffixMatcher.matches()) {
                Matcher shortSuffixMatcher = Pattern.compile("([a-zA-Z]*)(ene|yne|thiol|ol|imine|al|amine|phosphine|arsine|one|amide|oic acid|carboxylic acid|ether|oate|ane)").matcher(matcher.group(2));
                if (shortSuffixMatcher.matches()) {
                    pretoken = (radicalInMainChainFull = shortSuffixMatcher.group(1));
                    suftoken = shortSuffixMatcher.group(2);
                }
            } else {
                pretoken = suffixMatcher.group(1);
                suftoken = suffixMatcher.group(5);

                String firstRadicalOption = suffixMatcher.group(2);
                String secondRadicalOption = suffixMatcher.group(4);

                if (firstRadicalOption != null && !firstRadicalOption.isEmpty() || secondRadicalOption != null && !secondRadicalOption.isEmpty()) {
                    radicalInMainChainFull = firstRadicalOption != null && !firstRadicalOption.isEmpty() ? firstRadicalOption : secondRadicalOption;
                }

                if (secondRadicalOption != null && !secondRadicalOption.isEmpty()) {
                    Matcher ramificationMatcher = Pattern.compile("(?:(?:\\d+,?)+-)*[a-zA-Z]*(?:en|yn|enyn|an)").matcher(pretoken.substring(radicalInMainChainFull.length()));
                    while (ramificationMatcher.find()) {
                        String ramificationStr = ramificationMatcher.group(0);
                        if (!ramificationStr.endsWith("an")) {
                            posCount = ramificationStr.indexOf("-") < 0 ? 1 : (int) Arrays.stream(ramificationStr.substring(0, ramificationStr.indexOf("-")).split(",")).count();

                            if (posCount > 1) {
                                multiplier = MULTIPLIERS[posCount - 2];
                                if (!ramificationStr.substring(ramificationStr.indexOf("-") + 1).startsWith(multiplier)) {
                                    throw new RuntimeException("Multiplier of alkenes/alkynes doesn't match the number of position, molec: " + molec);
                                }
                            }

                            if (ramificationStr.endsWith("en")) {
                                numberOfAlkeneBounds += posCount;
                            } else if (ramificationStr.endsWith("yn") && !ramificationStr.endsWith("enyn")) {
                                numberOfAlkyneBounds += posCount;
                            } else if (ramificationStr.endsWith("enyn")) {
                                numberOfAlkeneBounds += posCount;
                                numberOfAlkyneBounds += posCount;
                            }
                        }
                    }
                } else {
                    String alkSuffix = suffixMatcher.group(3);
                    if (alkSuffix != null && !alkSuffix.isEmpty() && !alkSuffix.endsWith("an")) {
                        if (alkSuffix.endsWith("en")) {
                            numberOfAlkeneBounds += 1;
                        } else if (alkSuffix.endsWith("yn") && !alkSuffix.endsWith("enyn")) {
                            numberOfAlkyneBounds += 1;
                        } else if (alkSuffix.endsWith("enyn")) {
                            numberOfAlkeneBounds += 1;
                            numberOfAlkyneBounds += 1;
                        }
                    }
                }
                
                if ((pretoken == null || pretoken.isEmpty()) && (suftoken == null || suftoken.isEmpty())) {
                    suftoken = SUFFIXES[10]; // benzene
                }
            }

            if (pretoken != null && !pretoken.isEmpty()) {
                // The order of steps is a little strange. It may sound absurd when we process the subchain before the main chain, here I don't want to repeat the code and therefore change the order
                cyclic = pretoken.startsWith(PREFIXES[0]);
                if (cyclic) {
                    numberOfCyclos++;
                }
                String radicalInMainGroupAlmostFull = radicalInMainChainFull.substring(cyclic ? PREFIXES[0].length() : 0, radicalInMainChainFull.indexOf("-") < 0 ? radicalInMainChainFull.length() : radicalInMainChainFull.indexOf("-"));
                if (radicalInMainGroupAlmostFull.endsWith("an") || radicalInMainGroupAlmostFull.endsWith("en") || radicalInMainGroupAlmostFull.endsWith("yn")) {
                    if (radicalInMainGroupAlmostFull.endsWith("en")) {
                        numberOfAlkeneBounds += 1;
                    } else if (radicalInMainGroupAlmostFull.endsWith("yn") && !radicalInMainGroupAlmostFull.endsWith("enyn")) {
                        numberOfAlkyneBounds += 1;
                    } else if (radicalInMainGroupAlmostFull.endsWith("enyn")) {
                        numberOfAlkeneBounds += 1;
                        numberOfAlkyneBounds += 1;
                    }
                    radicalInMainGroupAlmostFull = radicalInMainGroupAlmostFull.substring(0, radicalInMainGroupAlmostFull.length() - 2);
                }
                final String radicalInMainGroup = radicalInMainGroupAlmostFull;
                optRadicalIdx = IntStream.range(0, RADICALS.length).boxed().filter(i -> radicalInMainGroup.equals(RADICALS[i])).findFirst();
        
                if (optRadicalIdx.isPresent()) {
                    numberOfCarbons += (optRadicalIdx.get() + 1);
                }
            }
            if (suftoken.endsWith("ane")) {
                // 4 * n - 2 * (n - 1) = 2 * n + 2 because each bound takes away 2 Hydrogens
                // Each cyclo takes 2 Hydrogens away
                result.put(CARBON, numberOfCarbons);
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith("ene") && !suftoken.endsWith(SUFFIXES[10])) { // not "benzene"
                posCount = suftoken.lastIndexOf("-") < 0 ? 1 : (int) Arrays.stream(suftoken.substring(1, suftoken.lastIndexOf("-")).split(",")).count();
                if (posCount > 1) {
                    multiplier = MULTIPLIERS[posCount - 2];

                    if (!suftoken.substring(suftoken.lastIndexOf("-") + 1, suftoken.indexOf("ene")).startsWith(multiplier)) {
                        throw new RuntimeException("Multiplier of alkenes doesn't match the number of position, molec: " + molec);
                    }
                }
                // 4 * n - 2 * (n - 1) = 2 * n + 2 because each bound takes away 2 Hydrogens
                // Each cyclo takes 2 Hydrogens away
                // Each Alkene bound takes away 2 Hydrogens
                // 2 * n + 2 - 2 * #cyclos - 2 * #AlkeneBounds
                result.put(CARBON, numberOfCarbons);
                numberOfAlkeneBounds += posCount;
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith("yne")) { // yne
                posCount = suftoken.lastIndexOf("-") < 0 ? 1 : (int) Arrays.stream(suftoken.substring(1, suftoken.lastIndexOf("-")).split(",")).count();
                if (posCount > 1) {
                    multiplier = MULTIPLIERS[posCount - 2];

                    if (!suftoken.substring(suftoken.lastIndexOf("-") + 1, suftoken.indexOf("yne")).startsWith(multiplier)) {
                        throw new RuntimeException("Multiplier of alkynes doesn't match the number of position, molec: " + molec);
                    }
                }
                // 4 * n - 2 * (n - 1) = 2 * n + 2 because each bound takes away 2 Hydrogens
                // Each cyclo takes 2 Hydrogens away
                // Each Alkene bound takes away 2 Hydrogens
                // Each Alkyne bound takes away 4 Hydrogens
                // 2 * n + 2 - 2 * #AlkeneBounds - 4 * #AlkyneBounds
                result.put(CARBON, numberOfCarbons);
                numberOfAlkyneBounds += posCount;
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith(SUFFIXES[11])) { // thiol
                posCount = parseSuffixTerm(suftoken, SUFFIXES[11]);
                // in case of SH: number of Hydrogens won't decrease by 1
                result.put("S", posCount);
                result.put(CARBON, numberOfCarbons);
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith(SUFFIXES[0])) { // ol
                posCount = parseSuffixTerm(suftoken, SUFFIXES[0]);
                // in case of OH: number of Hydrogens won't decrease by 1
                result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + posCount);
                result.put(CARBON, numberOfCarbons);
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith(SUFFIXES[9])) { // imine
                posCount = parseSuffixTerm(suftoken, SUFFIXES[9]);
                // in case of NH: each Carbon is losing 2 Hydrogens but instead gains 1 Hydro, so it'll lose 1 eventually
                result.put(NITROGEN, (result.get(NITROGEN) != null ? result.get(NITROGEN) : 0) + posCount);
                result.put(CARBON, numberOfCarbons);
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds - posCount);
            } else if (suftoken.endsWith(SUFFIXES[1]) || suftoken.endsWith(SUFFIXES[2])) { // al or one
                int alEnum = suftoken.endsWith(SUFFIXES[1]) ? 0 : 1;
                posCount = parseSuffixTerm(suftoken, alEnum == 0 ? SUFFIXES[1] : SUFFIXES[2]);
                // in case of =O: number of Hydrogens will decrease by 2 for each posCount
                result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + posCount);
                result.put(CARBON, numberOfCarbons);
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds - posCount) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith(SUFFIXES[8]) || suftoken.endsWith(SUFFIXES[12]) || suftoken.endsWith(SUFFIXES[13])) { // amine or phosphine or arsine
                int amineEnum = suftoken.endsWith(SUFFIXES[8]) ? 0 : (suftoken.endsWith(SUFFIXES[12]) ? 1 : 2);
                posCount = parseSuffixTerm(suftoken, amineEnum == 0 ? SUFFIXES[8] : (amineEnum == 1 ? SUFFIXES[12] : SUFFIXES[13]));
                String elem = "";
                switch (amineEnum) {
                    case 0:
                        elem = NITROGEN;
                        break;
                    case 1:
                        elem = "P";
                        break;
                    default:
                        elem = "As";                            
                }

                // in case of -NH2: each Carbon is losing one Hydro but gaining 2 Hydrogens, so #Hydrogens will increase by 1 for each posCount
                result.put(elem, (result.get(elem) != null ? result.get(elem) : 0) + posCount);
                result.put(CARBON, numberOfCarbons);
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds + posCount);
            } else if (suftoken.endsWith(SUFFIXES[7])) { // amide
                posCount = parseSuffixTerm(suftoken, SUFFIXES[7]);
                result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + posCount);
                result.put(NITROGEN, (result.get(NITROGEN) != null ? result.get(NITROGEN) : 0) + posCount);
                result.put(CARBON, numberOfCarbons);
                // each -NH2-C=O takes away 1 Hydrogen
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds - posCount);
            } else if (suftoken.endsWith(SUFFIXES[3]) || suftoken.endsWith(SUFFIXES[4])) { // "oic acid" or "carboxylic acid"
                int carboxylicEnum = suftoken.endsWith(SUFFIXES[3]) ? 0 : 1;
                posCount = parseSuffixTerm(suftoken, carboxylicEnum == 0 ? SUFFIXES[3] : SUFFIXES[4]);
                result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + 2 * posCount);
                numberOfCarbons += (carboxylicEnum == 0 ? 0 : 1);
                result.put(CARBON, numberOfCarbons);
                // each OH-C=O takes away 2 Hydrogens
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds - 2 * posCount);
            } else if (suftoken.endsWith(SUFFIXES[6])) { // "ether"
                posCount = parseSuffixTerm(suftoken, SUFFIXES[6]);
                result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + posCount);
                result.put(CARBON, numberOfCarbons);
                // each -O- doesn't change #Hydrogens
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith(SUFFIXES[5])) { // "oate"
                posCount = parseSuffixTerm(suftoken, SUFFIXES[5]);
                result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + 2 * posCount);
                result.put(CARBON, numberOfCarbons);
                // each oate function takes away 2 Hydrogens
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds - posCount) + 2 - numberOfFunctionBounds);
            } else if (suftoken.endsWith(SUFFIXES[10])) { // "benzene"
                posCount = parseSuffixTerm(suftoken, SUFFIXES[10]);
                // each benzene contains 6 Carbons
                numberOfCarbons += 6;
                result.put(CARBON, numberOfCarbons);
                // each benzene loses the same amount of Hydrogens as 3 alkene bounds and one Cyclo do altogether
                result.put(HYDROGEN, 2 * (numberOfCarbons - numberOfCyclos - numberOfAlkeneBounds - 2 * numberOfAlkyneBounds - 4 * posCount) + 2 - numberOfFunctionBounds);
            } else {
                throw new RuntimeException("Suffix " + suftoken + " has not been implemented, molec: " + molec);
            }

            result.entrySet().stream().filter(entry -> entry.getValue() == 0).map(entry -> entry.getKey()).collect(Collectors.toSet())
                    .forEach(key -> result.remove(key));

            return result;
        }

        return result;
    }

    private int parseSuffixTerm(String suftoken, String suffix) {
        int dashIdx = suftoken.indexOf("-");
        int posCount = 0;
        if (dashIdx > 0) {
            posCount = dashIdx < 0 ? 1 : (int) Arrays.stream(suftoken.substring(0, dashIdx).split(",")).count();
            if (posCount > 1) {
                String multiplier = MULTIPLIERS[posCount - 2];

                if (!suftoken.substring(dashIdx + 1).startsWith(multiplier)) {
                    throw new RuntimeException("Multiplier doesn't match the number of position, molec: " + molec);
                }                                    
            }
        } else {
            final String suftokenConstant = suftoken.substring(0, suftoken.indexOf(suffix));
            Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                .boxed()
                .sorted((i1, i2) -> i2 - i1)
                .filter(i -> suftokenConstant.equals(MULTIPLIERS[i]))
                .findFirst();

            posCount = optMultiplierIdx.orElse(-1) + 2;
        }

        return posCount;
    }

    private void parseSubchain(String token, int principalMultiple, RelativeIndexHolder spaceIdxHolder) {
        StringBuilder tokenBuilder = new StringBuilder(token);

        int openBrCounter = 1;
        int openBrIdx = token.indexOf("[");

        if (openBrIdx < 0) {
            boolean continueToParse = true;
            while (tokenBuilder.chars().boxed().filter(c -> c != ' ' && c != '-').count() > 0 && continueToParse) {
                if (tokenBuilder.charAt(0) == '-') {
                    tokenBuilder.deleteCharAt(0);
                }
                if (tokenBuilder.charAt(0) == ' ') {
                    spaceIdxHolder.setRelIdx(-1);
                    tokenBuilder.deleteCharAt(0);
                }
                continueToParse = parsePrefixStep(tokenBuilder, principalMultiple, spaceIdxHolder);
            }
            if (!continueToParse) {
                throw new RuntimeException(String.format("Parse exception while parsing subchain step, token: %s, molec: %s", tokenBuilder.toString(), molec));
            }
            return;
        }

        int lastDashIdx = token.substring(0, openBrIdx).lastIndexOf("-");
        Optional<Integer> optMultiplierIdx = null;
        int multiple = -1;
        int startIdx = lastDashIdx + 1;
        int spaceIdx = token.indexOf(" ");
        if (lastDashIdx >= 0) {
            int firstDashIdx = token.substring(0, lastDashIdx).lastIndexOf("-");
            int posCount = -1;
            if (firstDashIdx >= 0) {
                posCount = (int) Arrays.stream(token.substring(firstDashIdx + 1, lastDashIdx).split(",")).count();
            } else if (!Character.isDigit(token.charAt(0)) && !(spaceIdx >= 0 && spaceIdx < lastDashIdx)) {
                throw new RuntimeException(String.format("Dangle dash followed by subchain, molec: %s", molec));
            } else {
                posCount = (int) Arrays.stream(token.substring(Character.isDigit(token.charAt(0)) ? 0 : spaceIdx + 1, lastDashIdx).split(",")).count();
            }
            if (token.substring(lastDashIdx + 1, openBrIdx).equals(posCount > 1 ? MULTIPLIERS[posCount - 2] : "")) {
                multiple = posCount;
            } else if (posCount > 1 && !token.substring(lastDashIdx + 1).startsWith(MULTIPLIERS[posCount - 2])) {
                throw new RuntimeException("Multipliers of subchain doesn't match the number of proceeding positions, molec:" + molec);
            } else if (posCount > 1) {
                startIdx += MULTIPLIERS[posCount - 2].length();
            }
        }
        if (multiple < 0) {
            final String multiplierCand = token.substring(startIdx, openBrIdx);
            optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                .boxed()
                .sorted((i1, i2) -> i2 - i1)
                .filter(i -> multiplierCand.endsWith(MULTIPLIERS[i]))
                .findFirst();
        
            multiple = optMultiplierIdx.orElse(-1) + 2;
        }

        int charIdx = openBrIdx + 1;
        while (charIdx < token.length() && openBrCounter > 0) {
            if (token.charAt(charIdx) == '[') {
                ++openBrCounter;
            } else if (token.charAt(charIdx) == ']') {
                --openBrCounter;
            }
            charIdx++;
        }

        int closeBrIdx = charIdx - 1;
        String simpleStrippedToken = token.substring(openBrIdx + 1, closeBrIdx);
        RelativeIndexHolder newSpaceIdxHolder = new RelativeIndexHolder(spaceIdxHolder.getRelIdx());
        if (spaceIdxHolder.getRelIdx() >= 0) {
            newSpaceIdxHolder.updateInsideBrackets(openBrIdx, closeBrIdx);
        }
        if (!simpleStrippedToken.contains("[") && !simpleStrippedToken.contains("]")) {
            parsePrefix(simpleStrippedToken, multiple * principalMultiple, newSpaceIdxHolder);
        } else {
            parseSubchain(simpleStrippedToken, multiple * principalMultiple, newSpaceIdxHolder);
        }
        
        tokenBuilder.delete(openBrIdx, closeBrIdx + 1);

        if (tokenBuilder.isEmpty()) {
            return;
        }

        if (tokenBuilder.charAt(0) == '-') {
            tokenBuilder.deleteCharAt(0);
        }
        if (tokenBuilder.charAt(0) == ' ') {
            spaceIdxHolder.setRelIdx(-1);
            tokenBuilder.deleteCharAt(0);
        }
        newSpaceIdxHolder = new RelativeIndexHolder(spaceIdxHolder.getRelIdx());
        if (spaceIdxHolder.getRelIdx() >= 0) {
            newSpaceIdxHolder.updateOnShortenedString(openBrIdx, closeBrIdx);
        }
        parsePrefixStep(tokenBuilder, principalMultiple, newSpaceIdxHolder);

        if (!tokenBuilder.isEmpty()) {
            if (tokenBuilder.charAt(0) == '-') {
                tokenBuilder.deleteCharAt(0);
            }
            if (tokenBuilder.charAt(0) == ' ') {
                newSpaceIdxHolder.setRelIdx(-1);
                tokenBuilder.deleteCharAt(0);
            }
            if (!tokenBuilder.toString().contains("[") && !tokenBuilder.toString().contains("]")) {
                boolean continueToParse = true;
                while (tokenBuilder.chars().boxed().filter(c -> c != ' ' && c != '-').count() > 0 && continueToParse) {
                    if (tokenBuilder.charAt(0) == '-') {
                        tokenBuilder.deleteCharAt(0);
                    }
                    continueToParse = parsePrefixStep(tokenBuilder, principalMultiple, newSpaceIdxHolder);
                }
                if (!continueToParse) {
                    throw new RuntimeException(String.format("Parse exception while parsing subchain step, token: %s, molec: %s", tokenBuilder.toString(), molec));
                }
            } else {
                parseSubchain(tokenBuilder.toString(), principalMultiple, newSpaceIdxHolder);
            }
        }
    }

    private boolean parsePrefixStep(StringBuilder tokenBuilder, int multiple, RelativeIndexHolder spaceIdxHolder) {
        String currentToken = null;
        int openBrIdx = tokenBuilder.indexOf("[");
        String limitedToken = openBrIdx > 0 ? tokenBuilder.substring(0, openBrIdx) : tokenBuilder.toString();
        Matcher matcher = OYLOXY_PATTERN.matcher(limitedToken);
        String tmpToken = null;
        int tmpIdx = 0;

        // verify oyloxy syntax
        if (matcher.find() && validateOyloxy((tmpToken = matcher.group(0)), 0) 
                && tmpToken != null && tokenBuilder.indexOf(tmpToken) == 0 && ((tmpIdx = tokenBuilder.indexOf("[")) == -1 || tmpToken.length() <= tmpIdx)) {
            currentToken = tmpToken;
        } else {
            // !!! This is where your order is determined. The matcher.find() loop in parsePrefix() is therefore redundant
            // You've got to validate the syntaxes of both greedy and minium pattern 
            matcher = ALKYL_FUNCTION_PATTERN.matcher(limitedToken);
            // Prioritize ...oxy and check validity seems more feasible
            if (matcher.find() && validateAlkyl((tmpToken = matcher.group(0))) && tmpToken != null && ((tmpIdx = tokenBuilder.indexOf("[")) == -1 || tmpToken.length() <= tmpIdx)) {
                currentToken = tmpToken;
                Matcher oxyMatcher = OXY_PATTERN.matcher(tmpToken);
                if (oxyMatcher.find() && validateOyloxy(matcher.group(0), 1)) {
                    String oxyToken = oxyMatcher.group(0);
                    if (!currentToken.substring(currentToken.indexOf(oxyToken) + oxyToken.length()).startsWith("carbonyl")) {
                        currentToken = oxyToken;
                    }
                }
            } else {
                matcher = OXY_PATTERN.matcher(limitedToken);
                if (matcher.find() && validateOyloxy((tmpToken = matcher.group(0)), 1) && tmpToken != null && ((tmpIdx = tokenBuilder.indexOf("[")) == -1 || tmpToken.length() <= tmpIdx)) {
                    currentToken = tmpToken;
                }
            }
        }
        if (currentToken != null) {
            int oldSpaceIdx = spaceIdxHolder.getRelIdx();
            parsePrefix(currentToken, multiple, new RelativeIndexHolder(oldSpaceIdx)); // TODO: add [\\d]- when multiple == 1

            int startIdx = tokenBuilder.indexOf(currentToken);
            int endIdx = startIdx + currentToken.length();
            if (oldSpaceIdx >= 0) {
                spaceIdxHolder.updateOnShortenedString(startIdx, endIdx - 1);
            }
            tokenBuilder.delete(startIdx, endIdx);
            return true;
        }

        return false;
    }

    private final Pattern OYLOXY_PATTERN = Pattern.compile("(?:((?:\\d+,?)+-)?([a-zA-Z]*?)(en|yn|enyn|an)oyloxy|((?:\\d+,?)+-)?([a-zA-Z]*?)((?:(?:-((?:\\d+,?)+)-)?[a-zA-Z]*(?:en|yn|enyn))*)oyloxy)");

    private final Pattern ALKYL_FUNCTION_PATTERN = Pattern.compile("(?:(?:(?:\\d+,?)+-)?[a-zA-Z]*?(?:phenyl|formyl|hydroxy|bromo|imino|amino|phosphino|arsino|iodo|amido|mercapto|fluoro|chloro|oxo|carboxy)|([a-zA-Z]*?)yl|(((?:\\d+,?)+-)?([a-zA-Z]+?)((?:-(?:\\d+,?)+-[a-zA-Z]*(?:en|yn|enyn))*)(oxycarbon)?yl))-?");

    private final Pattern OXY_PATTERN = Pattern.compile("(?:((?:\\d+,?)+-)?([a-zA-Z]*?)(en|yn|enyn)?oxy|((?:\\d+,?)+-)?([a-zA-Z]*?)((?:(?:-((?:\\d+,?)+)-)?[a-zA-Z]*(?:en|yn|enyn))*)oxy)");

    private final Pattern PREFIX_PATTERN = Pattern.compile("((?:\\d+,?)+-)?[a-zA-Z]*(?:(?:(?:(?:-(?:\\d+,?)+-)?[a-zA-Z]*(?:en|yn|enyn))*|an)oyloxy|phenyl|formyl|(?:-(?:\\d+,?)+-[a-zA-Z]*(?:en|yn|enyn)|oxycarbon)*yl|mo|imino|amino|phosphino|arsino|amido|iodo|to|fluoro|chloro|oxo|carboxy|hydroxy|(?:(?:(?:-(?:\\d+,?)+-)?[a-zA-Z]*(?:en|yn|enyn))*)oxy)(?:-| )?");

    private boolean validateOyloxy(String token, int oyloxyEnum) {
        Matcher matcher = (oyloxyEnum == 0 ? OYLOXY_PATTERN : OXY_PATTERN).matcher(token);

        int multiple = 1;
        if (matcher.find()) {
            String currentWholeWord = matcher.group(0);

            if (matcher.group(2) != null && !matcher.group(2).isEmpty()) {
                int radical = 0;
                String pretoken = matcher.group(1);
                final String shortAlk = matcher.group(2);

                if (pretoken == null || pretoken.isEmpty()) {
                    // !!! possible ambiguous case: tridecoxy
                    Optional<Integer> optRadicalIdx = IntStream.range(0, RADICALS.length)
                        .boxed()
                        .sorted((i1, i2) -> i2 - i1)
                        .filter(i -> shortAlk.endsWith(RADICALS[i]))
                        .findFirst();
                    if (optRadicalIdx.isPresent()) {
                        if ("eth".equals(RADICALS[optRadicalIdx.get()]) && shortAlk.contains("meth")) {
                            optRadicalIdx = Optional.of(0);
                        }
                        radical = optRadicalIdx.get() + 1;
                        int lastIdx = shortAlk.length() - RADICALS[optRadicalIdx.get()].length();
                        boolean cyclic = shortAlk.substring(0, lastIdx).endsWith(PREFIXES[0]);
                        if (cyclic) {
                            lastIdx -= PREFIXES[0].length();
                        }
                        final String multiplierCand = shortAlk.substring(0, lastIdx);
                        Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                            .boxed()
                            .filter(i -> multiplierCand.equals(MULTIPLIERS[i]))
                            .findFirst();
                        if (optMultiplierIdx.isEmpty() && !multiplierCand.isEmpty()) {
                            return false;
                        }
                        multiple = optMultiplierIdx.orElse(-1) + 2;
                        return shortAlk.equals((multiple > 1 ? MULTIPLIERS[multiple - 2] : "") + (cyclic ? PREFIXES[0] : "") + RADICALS[radical - 1]);
                    }
                    return false;
                } else {
                    multiple = (int) Arrays.stream(pretoken.substring(0, pretoken.lastIndexOf("-")).split(",")).count();
                    if (multiple > 1) {
                        if (!shortAlk.startsWith(MULTIPLIERS[multiple - 2])) {
                            return false;
                        }
                    }
                    int lastIdx = multiple > 1 ? MULTIPLIERS[multiple - 2].length() : 0;
                    boolean cyclic = shortAlk.substring(lastIdx).startsWith(PREFIXES[0]);
                    final int lastIdxConst = cyclic ? lastIdx + PREFIXES[0].length() : lastIdx;
                    Optional<Integer> optRadicalIdx = IntStream.range(0, RADICALS.length)
                        .boxed()
                        .sorted((i1, i2) -> i2 - i1)
                        .filter(i -> shortAlk.substring(lastIdxConst).equals(RADICALS[i]))
                        .findFirst();

                    return optRadicalIdx.isPresent();
                }
            } else {
                if (matcher.group(4) == null && matcher.group(5) == null && matcher.group(6) == null) {
                    return false;
                }
                String pretoken = matcher.group(4);
                int posCount = (int) Arrays.stream(pretoken.substring(0, pretoken.lastIndexOf("-")).split(",")).count();

                String inBetweenWord = matcher.group(5);
                int startIdx = 0;
                String multiplier = null;
                if (posCount > 1) {
                    multiplier = MULTIPLIERS[posCount - 2];
                    startIdx = multiplier.length();
                    if (!inBetweenWord.startsWith(multiplier)) {
                        return false;
                    }
                }
                boolean cyclic = inBetweenWord.substring(posCount > 1 ? multiplier.length() : 0).startsWith(PREFIXES[0]);
                if (cyclic) {
                    startIdx += PREFIXES[0].length();
                }
                final String radicalCand = inBetweenWord.substring(startIdx);
                Optional<Integer> optRadicalIdx = IntStream.range(0, RADICALS.length)
                    .boxed()
                    .filter(i -> radicalCand.equals(RADICALS[i]))
                    .findFirst();

                multiple = posCount;
                String suftoken = matcher.group(6);
                if (optRadicalIdx.isPresent() && suftoken != null && !suftoken.isEmpty()) {
                    return validateRamification(suftoken);
                }
                return false;
            }
        }
        return false;
    }

    private boolean validateRamification(String ramificationStr) {
        Matcher alkMatcher = Pattern.compile("-(?:\\d+,?)+-[a-zA-Z]*(?:en|yn|enyn)").matcher(ramificationStr);
        while (alkMatcher.find()) {
            String suftoken = alkMatcher.group(0);
            int posCount = (int) Arrays.stream(suftoken.substring(1, suftoken.lastIndexOf("-")).split(",")).count();
            int startIdx = suftoken.lastIndexOf("-") + 1;
            if (posCount > 1) {
                String multiplier = MULTIPLIERS[posCount - 2];
                //startIdx = suftoken.lastIndexOf("-") + 1;
                if (!suftoken.substring(startIdx).startsWith(multiplier)) {
                    return false;
                }

                startIdx += multiplier.length();
            }
            String lastToken = suftoken.substring(startIdx);
            if (!"en".equals(lastToken) && !"yn".equals(lastToken)) {
                return false;
            }
        }
        return true;
    }

    private boolean validateAlkyl(String token) {
        if ((token.endsWith("yl") || token.endsWith("yl-")) && !token.endsWith("formyl") && !token.endsWith("formyl-") 
                && !token.endsWith("phenyl") && !token.endsWith("phenyl-")) {
            return alkylPattern.matcher(token).matches();
        }

        if (functionPattern.matcher(token).matches()) {
            return true;
        }
        Optional<String> optPrefix = Arrays.stream(PREFIXES)
            .filter(prefix -> (token.endsWith(prefix) || token.endsWith(prefix + "-")) && !PREFIXES[0].equals(prefix))
            .findFirst();

        if (optPrefix.isPresent()) {
            StringBuilder tokenBuilder = new StringBuilder(token.substring(0, token.indexOf(optPrefix.get())));
            String shortenedToken = tokenBuilder.toString();
            
            Matcher matcher;
            if ((matcher = alkylPattern.matcher(token)).find()) {
                String matchedAlkyl = matcher.group(0);
                String multiplierCand = token.substring(token.indexOf(matchedAlkyl) + matchedAlkyl.length(), token.indexOf(optPrefix.get()));
                Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                    .boxed()
                    .filter(i -> multiplierCand.equals(MULTIPLIERS[i]))
                    .findFirst();
                return optMultiplierIdx.isPresent();
            }
            return validateOyloxy(shortenedToken, 0) || validateAlkyl(shortenedToken) || validateOyloxy(shortenedToken, 1);
        }
        return false;
    }

    private void parsePrefix(String token, int principalMultiple, RelativeIndexHolder spaceIdxHolder) {
        int posCount = 0;
        String multiplier = null;
        Optional<Integer> optRadicalIdx;
        int startIdx = 0;
        int lastIdx = 0;
        boolean cyclic = false;
        String matchedWholeWord = null;

        Matcher prefixMatcher = PREFIX_PATTERN.matcher(token);

        while (prefixMatcher.find()) {
            // consume currentWholeWord
            Matcher oyloxyMatcher = OYLOXY_PATTERN.matcher(prefixMatcher.group(0));
            int multiple = 1;
            if (oyloxyMatcher.find()) {
                parsePrefixOyloxy(oyloxyMatcher, principalMultiple, spaceIdxHolder, 0);
            } else {
                StringBuilder ramificationBuilder = new StringBuilder(prefixMatcher.group(0));
                Matcher ramificationMatcher = ALKYL_FUNCTION_PATTERN.matcher(ramificationBuilder.toString());
                while (ramificationMatcher.find()) {
                    matchedWholeWord = ramificationMatcher.group(0);
                    int spaceIdx = spaceIdxHolder.getRelIdx();
                    boolean esterSuffix = esterSuffixMultiplier > 1 && spaceIdx >= matchedWholeWord.length(); // "oate"

                    if ((matchedWholeWord.endsWith("yl") || matchedWholeWord.endsWith("yl-")) && !matchedWholeWord.endsWith("formyl") && !matchedWholeWord.endsWith("formyl-") 
                            && !matchedWholeWord.endsWith("phenyl") && !matchedWholeWord.endsWith("phenyl-")) {
                        if (ramificationMatcher.group(2) != null && !ramificationMatcher.group(2).isEmpty()) {
                            String pretoken = ramificationMatcher.group(3);
                            String inBetweenToken = ramificationMatcher.group(4);
                            startIdx = 0;
                            if (pretoken != null && !pretoken.isEmpty()) {
                                posCount = (int) Arrays.stream(pretoken.substring(0, pretoken.indexOf("-")).split(",")).count();
                
                                if (posCount > 1) {
                                    //pos = Integer.parseInt(position);
                                    multiplier = MULTIPLIERS[posCount - 2];

                                    if (!inBetweenToken.startsWith(multiplier)) {
                                        throw new RuntimeException("Multiplier in ALKYL_PATTERN doesn't match the number  position, molec: " + molec);
                                    }
                                }
                                if (posCount > 1) {
                                    startIdx = MULTIPLIERS[posCount - 2].length();
                                }
                                cyclic = inBetweenToken.substring(startIdx).startsWith(PREFIXES[0]);
                                if (cyclic) {
                                    startIdx += PREFIXES[0].length();
                                }
                                lastIdx = inBetweenToken.length();
                                if (inBetweenToken.endsWith("en")) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkeneBounds += posCount * principalMultiple;
                                    } else {
                                        numberOfAlkeneBounds += posCount * principalMultiple;
                                    }
                                    lastIdx = inBetweenToken.lastIndexOf("en");
                                } else if (inBetweenToken.endsWith("yn") && !inBetweenToken.endsWith("enyn")) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkyneBounds += posCount * principalMultiple;
                                    } else {
                                        numberOfAlkyneBounds += posCount * principalMultiple;
                                    }
                                    lastIdx = inBetweenToken.lastIndexOf("yn");
                                } else if (inBetweenToken.endsWith("enyn")) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkeneBounds += posCount * principalMultiple;
                                        esterSuffixNumberOfAlkyneBounds += posCount;
                                    } else {
                                        numberOfAlkeneBounds += posCount * principalMultiple;
                                        numberOfAlkyneBounds += posCount;
                                    }
                                    lastIdx = inBetweenToken.lastIndexOf("enyn");
                                }
                                final String radicalInRamification = inBetweenToken.substring(startIdx, lastIdx);

                                optRadicalIdx = IntStream.range(0, RADICALS.length)
                                    .boxed()
                                    .filter(i -> radicalInRamification.equals(RADICALS[i]))
                                    .findFirst();
                            } else {
                                // difficult case where multiplier is follwed by radical without the preceding digit sequence, for example dioctadec-12-enyl
                                final String multiplierAndRadical = inBetweenToken.endsWith("an") || inBetweenToken.endsWith("en") || inBetweenToken.endsWith("yn") && !inBetweenToken.endsWith("enyn")
                                    ? inBetweenToken.substring(0, inBetweenToken.length() - 2) 
                                    : inBetweenToken.endsWith("enyn") ? inBetweenToken.substring(inBetweenToken.length() - 4) : inBetweenToken ;
                                optRadicalIdx = IntStream.range(0, RADICALS.length)
                                    .boxed()
                                    .sorted((i1, i2) -> i2 - i1)
                                    .filter(i -> multiplierAndRadical.endsWith(RADICALS[i]))
                                    .findFirst();

                                if (optRadicalIdx.isPresent() && "eth".equals(RADICALS[optRadicalIdx.get()]) && multiplierAndRadical.endsWith("meth")) {
                                    optRadicalIdx = Optional.of(0);
                                }

                                final int radicalIdx = optRadicalIdx.get();
                                lastIdx = multiplierAndRadical.indexOf(RADICALS[radicalIdx]);
                                cyclic = multiplierAndRadical.substring(0, lastIdx).endsWith(PREFIXES[0]);
                                if (cyclic) {
                                    lastIdx -= PREFIXES[0].length();
                                }
                                final String multiplierCand = multiplierAndRadical.substring(0, lastIdx);
                                Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                                        .boxed()
                                        .sorted((i1, i2) -> i2 - i1)
                                        .filter(i -> multiplierCand.equals(MULTIPLIERS[i]))
                                        .findFirst();
                                posCount = optMultiplierIdx.orElse(-1) + 2;
                                if (inBetweenToken.endsWith("en")) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkeneBounds += posCount * principalMultiple;
                                    } else {
                                        numberOfAlkeneBounds += posCount * principalMultiple;
                                    }
                                } else if (inBetweenToken.endsWith("yn") && !inBetweenToken.endsWith("enyn")) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkyneBounds += posCount * principalMultiple;
                                    } else {
                                        numberOfAlkyneBounds += posCount * principalMultiple;
                                    }
                                } else if (inBetweenToken.endsWith("enyn")) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkeneBounds += posCount * principalMultiple;
                                        esterSuffixNumberOfAlkyneBounds += posCount;
                                    } else {
                                        numberOfAlkeneBounds += posCount * principalMultiple;
                                        numberOfAlkyneBounds += posCount;
                                    }
                                }
                            }

                            if (optRadicalIdx.isPresent()) {
                                if (esterSuffix) {
                                    esterSuffixNumberOfCarbons += (optRadicalIdx.get() + 1) * posCount * principalMultiple;
                                } else {
                                   // alkyl group doesn't effect the number of Hydrogens
                                    numberOfCarbons += (optRadicalIdx.get() + 1) * posCount * principalMultiple;
                                }
                                if (cyclic) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfCyclos += posCount * principalMultiple;
                                    } else {
                                        numberOfCyclos += posCount * principalMultiple;
                                    }
                                }
                                String suftoken = ramificationMatcher.group(5);
                                multiple = posCount;
                                if (suftoken != null && !suftoken.isEmpty()) {
                                    parseRamification(suftoken, multiple * principalMultiple, esterSuffix);
                                }
                                String oxycarbonylStr = ramificationMatcher.group(6); 
                                if (oxycarbonylStr != null && !oxycarbonylStr.isEmpty()) {
                                    numberOfCarbons += multiple * principalMultiple;
                                    result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + 2 * multiple * principalMultiple);
                                    numberOfAlkeneBounds += multiple * principalMultiple; // one oxycarbon- takes away 2 Hydrogens                                    
                                }
                            }
                        } else if (ramificationMatcher.group(1) != null && !ramificationMatcher.group(1).isEmpty()) { // special case for amines, phosphines, arsines in the form of ((multiplier)*(radical)yl)*amine for example
                            StringBuilder shortRamificationBuilder = new StringBuilder(ramificationMatcher.group(1));
                            int lastIdxOfOxycarbon = shortRamificationBuilder.lastIndexOf("oxycarbon");
                            boolean isOxycarbon = shortRamificationBuilder.toString().endsWith("oxycarbon");
                            if (isOxycarbon) {
                                shortRamificationBuilder.delete(lastIdxOfOxycarbon, shortRamificationBuilder.length());
                            }

                            boolean isAlkene = shortRamificationBuilder.toString().endsWith("en");
                            boolean isAlkyne = shortRamificationBuilder.toString().endsWith("yn") && !shortRamificationBuilder.toString().endsWith("enyn");
                            boolean isAlkenyne = shortRamificationBuilder.toString().endsWith("enyn");
                            if (isAlkene) {
                                shortRamificationBuilder.delete(shortRamificationBuilder.length() - 2, shortRamificationBuilder.length());                
                            } else {
                                if (isAlkyne) {
                                    shortRamificationBuilder.delete(shortRamificationBuilder.length() - 2, shortRamificationBuilder.length());
                                } else if (isAlkenyne) {
                                    shortRamificationBuilder.delete(shortRamificationBuilder.length() - 4, shortRamificationBuilder.length());
                                }
                            }

                            final String shortAlkyl = shortRamificationBuilder.toString();

                            optRadicalIdx = IntStream.range(0, RADICALS.length)
                                .boxed()
                                .sorted((i1, i2) -> i2 - i1)
                                .filter(i -> shortAlkyl.endsWith(RADICALS[i]))
                                .findFirst();

                            multiple = 1;
                            if (optRadicalIdx.isPresent()) {
                                if ("eth".equals(RADICALS[optRadicalIdx.get()]) && shortAlkyl.contains("meth")) {
                                    optRadicalIdx = Optional.of(0);
                                }

                                startIdx = shortAlkyl.length() - RADICALS[optRadicalIdx.get()].length();
                                if (startIdx > 0) {
                                    cyclic = shortAlkyl
                                        .substring(0, startIdx)
                                            .endsWith(PREFIXES[0]);
                                    if (cyclic) {
                                        startIdx -= PREFIXES[0].length();
                                    }

                                    if (startIdx > 0) {
                                        final String multiplierCand = shortAlkyl.substring(0, startIdx);
                                        Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                                            .boxed()
                                            .sorted((i1, i2) -> i2 - i1)
                                            .filter(i -> multiplierCand.endsWith(MULTIPLIERS[i]))
                                            .findFirst();
                                        if (optMultiplierIdx.isPresent()) {
                                            multiple = optMultiplierIdx.get() + 2;
                                        }
                                    }
                                }
                                if (cyclic) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfCyclos += multiple * principalMultiple;
                                    } else {
                                        numberOfCyclos += multiple * principalMultiple;
                                    }
                                }
                                int n = (isOxycarbon ? optRadicalIdx.get() + 2 : optRadicalIdx.get() + 1) * multiple * principalMultiple;
                                if (esterSuffix) {
                                    esterSuffixNumberOfCarbons += n;
                                } else {
                                    numberOfCarbons += n;
                                }
                                if (isOxycarbon) {
                                    if (esterSuffix) {
                                        esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + 2 * multiple * principalMultiple);
                                        // each O=C-O takes away the same number of Hydrogens as one alkene group does
                                        esterSuffixNumberOfAlkeneBounds += multiple * principalMultiple;                                        
                                    } else {
                                        result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + 2 * multiple * principalMultiple);
                                        // each O=C-O takes away the same number of Hydrogens as one alkene group does
                                        numberOfAlkeneBounds += multiple * principalMultiple;
                                    }
                                }
                                if (isAlkene) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkeneBounds += multiple * principalMultiple;
                                    } else {
                                        numberOfAlkeneBounds += multiple * principalMultiple;
                                    }
                                } else if (isAlkyne) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkyneBounds += multiple * principalMultiple;
                                    } else {
                                        numberOfAlkyneBounds += multiple * principalMultiple;
                                    }
                                } else if (isAlkenyne) {
                                    if (esterSuffix) {
                                        esterSuffixNumberOfAlkeneBounds += multiple * principalMultiple;
                                        esterSuffixNumberOfAlkyneBounds += multiple;
                                    } else {
                                        numberOfAlkeneBounds += multiple * principalMultiple;
                                        numberOfAlkyneBounds += multiple;
                                    }
                                }
                            }
                        }
                    } else {
                        Matcher prefixFunctionMatcher = Pattern.compile("((?:(?:\\d+,?)+-)?)([a-zA-Z]*(?:yl|mo|no|do|to|fluoro|chloro|bromo|iodo|xo|carboxy|hydroxy))-?").matcher(matchedWholeWord);
                        if (prefixFunctionMatcher.find()) {
                            String pretoken = prefixFunctionMatcher.group(1);
                            String suftoken = prefixFunctionMatcher.group(2);
                            String prefixFunctionFullString = prefixFunctionMatcher.group(0);

                            posCount = -1;
                            Optional<String> optPrefix = Arrays.stream(PREFIXES)
                                    .filter(prefix -> suftoken.endsWith(prefix) && !PREFIXES[0].equals(prefix))
                                    .findFirst();
                            
                            if (optPrefix.isPresent()) {
                                boolean suftokenPure = functionPattern.matcher(prefixFunctionFullString).matches();
                                if (suftokenPure) {
                                    if (pretoken != null && !pretoken.isEmpty()) {
                                        posCount = (int) Arrays.stream(pretoken.split(",")).count();

                                        if (posCount > 1) {
                                            multiplier = MULTIPLIERS[posCount - 2];

                                            if (!suftoken.startsWith(multiplier)) {
                                                throw new RuntimeException("Multiplier in PREFIX_FUNCTION_PATTERN doesn't match the number of position, molec: " + molec);
                                            }
                                        }
                                    } else {
                                        final String multiplierCand = suftoken.substring(0, suftoken.indexOf(optPrefix.get()));
                                        if ("".equals(multiplierCand)) {
                                            posCount = 1;
                                        } else {
                                            Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                                                .boxed()
                                                .filter(i -> multiplierCand.equals(MULTIPLIERS[i]))
                                                .findFirst();

                                            if (optMultiplierIdx.isPresent()) {
                                                posCount = optMultiplierIdx.get() + 2;
                                            }
                                        }  
                                    }
                                }

                                // example of the second predicate: 3,10,15-trihexoxy-4,6-diphenyl, 9-pent-2,3-diynylphenyl-
                                if (!suftokenPure || prefixFunctionFullString.length() < matchedWholeWord.length()) {
                                    int shortFunctionIdx = matchedWholeWord.indexOf(!suftokenPure ? optPrefix.get() : prefixFunctionFullString);
                                    String prefixFunctionShort = matchedWholeWord.substring(0, shortFunctionIdx);
                                    RelativeIndexHolder newSpaceIdxHolder = new RelativeIndexHolder(spaceIdxHolder.getRelIdx());
                                    newSpaceIdxHolder.updateOnShortenedString(shortFunctionIdx, shortFunctionIdx + (!suftokenPure ? optPrefix.get() : prefixFunctionFullString).length() - 1);
                                    // TODO: 9-heptyldihydroxy-
                                    StringBuilder tokenBuilder = new StringBuilder(prefixFunctionShort);
                                    boolean continueToParse = true;
                                    //parsePrefix(prefixFunctionFullString.substring(0, shortFunctionIdx), principalMultiple, newSpaceIdxHolder);
                                    while (continueToParse && tokenBuilder.chars().boxed().filter(c -> c != ' ' && c != '-').count() > 0) {
                                        if (tokenBuilder.charAt(0) == '-') {
                                                tokenBuilder.deleteCharAt(0);
                                        }
                                        if (tokenBuilder.charAt(0) == ' ') {
                                            spaceIdxHolder.setRelIdx(-1);
                                            tokenBuilder.deleteCharAt(0);
                                        }
                                        continueToParse = parsePrefixStep(tokenBuilder, principalMultiple, newSpaceIdxHolder);
                                    }
                                }
                                parseFunctionPrefix(optPrefix.get(), (suftokenPure || (!prefixFunctionFullString.equals(matchedWholeWord) && posCount > 0) ? posCount : 1) * principalMultiple, esterSuffix);
                            }
                        }
                    }
                    int ramificationStartIdx = ramificationBuilder.indexOf(matchedWholeWord);
                    int ramificationEndIdx = ramificationStartIdx + matchedWholeWord.length();
                    spaceIdxHolder.updateOnShortenedString(ramificationStartIdx, ramificationEndIdx - 1);
                    ramificationBuilder.delete(ramificationStartIdx, ramificationEndIdx);
                }

                if (!ramificationBuilder.isEmpty()) {
                    Matcher oxyMatcher = OXY_PATTERN.matcher(ramificationBuilder);
                    String currentWholeWord = null;
                    if (oxyMatcher.find() && !(currentWholeWord = oxyMatcher.group(0)).endsWith(PREFIXES[1]) && !currentWholeWord.endsWith(PREFIXES[3])) {
                        parsePrefixOyloxy(oxyMatcher, principalMultiple, spaceIdxHolder, 1);
                    }
                }
            } // end if 
        } // end while 
    }

    private void parseRamification(String ramificationStr, int multiple, boolean esterSuffix) {
        Matcher alkMatcher = Pattern.compile("-(?:\\d+,?)+-[a-zA-Z]*(?:en|yn|enyn)").matcher(ramificationStr);
        while (alkMatcher.find()) {
            String suftoken = alkMatcher.group(0);
            int posCount = (int) Arrays.stream(suftoken.substring(1, suftoken.lastIndexOf("-")).split(",")).count();
            int startIdx = suftoken.lastIndexOf("-") + 1;
            if (posCount > 1) {
                String multiplier = MULTIPLIERS[posCount - 2];
                //startIdx = suftoken.lastIndexOf("-") + 1;
                if (!suftoken.substring(startIdx).startsWith(multiplier)) {
                    throw new RuntimeException("Multiplier in ALKYL_PATTERN doesn't match the number  position, molec: " + molec);
                }

                startIdx += multiplier.length();
            }
            if ("en".equals(suftoken.substring(startIdx))) {
                if (esterSuffix) {
                    esterSuffixNumberOfAlkeneBounds += multiple * posCount;
                } else {
                    numberOfAlkeneBounds += multiple * posCount;
                }
            } else if ("yn".equals(suftoken.substring(startIdx)) && !"enyn".equals(suftoken.substring(startIdx))) {
                if (esterSuffix) {
                    esterSuffixNumberOfAlkyneBounds += multiple * posCount;
                } else {
                    numberOfAlkyneBounds += multiple * posCount;
                }
            } else if ("enyn".equals(suftoken.substring(startIdx))) {
                if (esterSuffix) {
                    esterSuffixNumberOfAlkeneBounds += multiple * posCount;
                    esterSuffixNumberOfAlkyneBounds += multiple;
                } else {
                    numberOfAlkeneBounds += multiple * posCount;
                    numberOfAlkyneBounds += multiple;
                }
            }
        }
    }

    private void parsePrefixOyloxy(Matcher matcher, int principalMultiple, RelativeIndexHolder spaceIdxHolder, int oyloxyEnum) {
        int multiple = 1;
        String currentWholeWord = matcher.group(0);
        int spaceIdx = spaceIdxHolder.getRelIdx();
        boolean esterSuffix = esterSuffixMultiplier > 1 && spaceIdx >= currentWholeWord.length(); // "oate" suffix

        if (matcher.group(2) != null && !matcher.group(2).isEmpty()) {
            int radical = 0;
            boolean cyclic = false;
            Optional<Integer> optRadicalIdx = null;
            String pretoken = matcher.group(1);
            final String shortAlk = matcher.group(2);

            if (pretoken == null || pretoken.isEmpty()) {
                optRadicalIdx = IntStream.range(0, RADICALS.length)
                    .boxed()
                    .sorted((i1, i2) -> i2 - i1)
                    .filter(i -> shortAlk.endsWith(RADICALS[i]))
                    .findFirst();

                if (optRadicalIdx.isPresent()) {
                    if ("eth".equals(RADICALS[optRadicalIdx.get()]) && shortAlk.contains("meth")) {
                        optRadicalIdx = Optional.of(0);
                    }
                    radical = optRadicalIdx.get() + 1;
                    int lastIdx = shortAlk.length() - RADICALS[optRadicalIdx.get()].length();
                    cyclic = shortAlk.substring(0, lastIdx).endsWith(PREFIXES[0]);
                    if (cyclic) {
                        lastIdx -= PREFIXES[0].length();
                    }
                    final String multiplierCand = shortAlk.substring(0, lastIdx);
                    Optional<Integer> optMultiplierIdx = IntStream.range(0, MULTIPLIERS.length)
                        .boxed()
                        .filter(i -> multiplierCand.equals(MULTIPLIERS[i]))
                        .findFirst();
                    multiple = optMultiplierIdx.orElse(-1) + 2;
                }
            } else {
                multiple = (int) Arrays.stream(pretoken.substring(0, pretoken.lastIndexOf("-")).split(",")).count();
                if (multiple > 1) {
                    if (!shortAlk.startsWith(MULTIPLIERS[multiple - 2])) {
                        throw new RuntimeException("Multiplier in OYLOXY_PATTERN doesn't match the number of position, molec: " + molec);
                    }
                }
                int lastIdx = multiple > 1 ? MULTIPLIERS[multiple - 2].length() : 0;
                cyclic = shortAlk.substring(lastIdx).startsWith(PREFIXES[0]);
                if (cyclic) {
                    lastIdx += PREFIXES[0].length();
                }
                final String radicalCand = shortAlk.substring(lastIdx);
                optRadicalIdx = IntStream.range(0, RADICALS.length)
                    .boxed()
                    .sorted((i1, i2) -> i2 - i1)
                    .filter(i -> radicalCand.equals(RADICALS[i]))
                    .findFirst();
            }
            if (optRadicalIdx.isPresent()) {
                if ("eth".equals(RADICALS[optRadicalIdx.get()]) && shortAlk.contains("meth")) {
                    radical = 1;
                } else {
                    radical = optRadicalIdx.get() + 1;
                }
                if (esterSuffix) {
                    esterSuffixNumberOfCarbons += radical * multiple * principalMultiple;
                } else {
                    numberOfCarbons += radical * multiple * principalMultiple;
                }
                if (cyclic) {
                    if (esterSuffix) {
                        esterSuffixNumberOfCyclos += multiple * principalMultiple;
                    } else {
                        numberOfCyclos += multiple * principalMultiple;
                    }
                }
                if ("en".equals(matcher.group(3))) {
                    if (esterSuffix) {
                        esterSuffixNumberOfAlkeneBounds += multiple * principalMultiple;
                    } else {
                        numberOfAlkeneBounds += multiple * principalMultiple;
                    }
            
                } else if ("yn".equals(matcher.group(3)) && !"enyn".equals(matcher.group(3))) {
                    if (esterSuffix) {
                        esterSuffixNumberOfAlkyneBounds += multiple * principalMultiple;
                    } else {
                        numberOfAlkyneBounds += multiple * principalMultiple;
                    }
                } else if ("enyn".equals(matcher.group(3))) {
                    if (esterSuffix) {
                        esterSuffixNumberOfAlkeneBounds += multiple * principalMultiple;
                        esterSuffixNumberOfAlkyneBounds += multiple;
                    } else {
                        numberOfAlkeneBounds += multiple * principalMultiple;
                        numberOfAlkyneBounds += multiple;
                    }
                }
            }
        } else {
            String pretoken = matcher.group(4);
            int posCount = (int) Arrays.stream(pretoken.substring(0, pretoken.lastIndexOf("-")).split(",")).count();
            int startIdx = 0;
            String inBetweenWord = matcher.group(5);
            if (posCount > 1) {
                String multiplier = MULTIPLIERS[posCount - 2];
                startIdx = multiplier.length();
                if (!inBetweenWord.startsWith(multiplier)) {
                    throw new RuntimeException("Multiplier in OYLOXY_PATTERN doesn't match the number of position, molec: " + molec);
                }
            }
            boolean cyclic = inBetweenWord.substring(startIdx).startsWith(PREFIXES[0]);
            if (cyclic) {
                startIdx += PREFIXES[0].length();
                if (esterSuffix) {
                    esterSuffixNumberOfCyclos += posCount * principalMultiple;
                } else {
                    numberOfCyclos += posCount * principalMultiple;
                }
            }
            final String radicalCand = inBetweenWord.substring(startIdx);
            Optional<Integer> optRadicalIdx = IntStream.range(0, RADICALS.length)
                .boxed()
                .filter(i -> radicalCand.equals(RADICALS[i]))
                .findFirst();

            multiple = posCount;
            String suftoken = matcher.group(6);
            if (optRadicalIdx.isPresent() && suftoken != null && !suftoken.isEmpty()) {
                int n = (optRadicalIdx.get() + 1) * multiple * principalMultiple;
                if (esterSuffix) {
                    esterSuffixNumberOfCarbons += n;
                } else {
                    numberOfCarbons += n;
                }
                parseRamification(suftoken, multiple * principalMultiple, esterSuffix);
            }
        }
        if (esterSuffix) {
            esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + (oyloxyEnum == 0 ? 2 : 1) * multiple * principalMultiple);
        } else {
            result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + (oyloxyEnum == 0 ? 2 : 1) * multiple * principalMultiple);
        }
        // alkanoyloxy in general takes away 2 Hydrogens
        // whereas alkoxy doesn't take away any Hydrogens
        if (esterSuffix && oyloxyEnum == 0) {
            esterSuffixNumberOfAlkeneBounds += multiple * principalMultiple;
        } else if (oyloxyEnum == 0) {
            numberOfAlkeneBounds += multiple * principalMultiple;
        }
    }

    private void parseFunctionPrefix(String suftoken, int multiple, boolean esterSuffix) {
        switch (suftoken) {
            case "fluoro":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put("F", (esterSuffixResult.get("F") != null ? esterSuffixResult.get("F") : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds += multiple;
                } else {
                    result.put("F", (result.get("F") != null ? result.get("F") : 0) + multiple);
                    numberOfFunctionBounds += multiple;
                }
                break;
            case "chloro":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put("Cl", (esterSuffixResult.get("Cl") != null ? esterSuffixResult.get("Cl") : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds += multiple;
                } else {
                    result.put("Cl", (result.get("Cl") != null ? result.get("Cl") : 0) + multiple);
                    numberOfFunctionBounds += multiple;
                }
                break;
            case "bromo":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put("Br", (esterSuffixResult.get("Br") != null ? esterSuffixResult.get("Br") : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds += multiple;
                } else {                
                    result.put("Br", (result.get("Br") != null ? result.get("Br") : 0) + multiple);
                    numberOfFunctionBounds += multiple;
                }
                break;
            case "iodo":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put("I", (esterSuffixResult.get("I") != null ? esterSuffixResult.get("I") : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds += multiple;
                } else {
                    result.put("I", (result.get("I") != null ? result.get("I") : 0) + multiple);
                    numberOfFunctionBounds += multiple;
                }
                break;
            case "hydroxy":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    // in case of OH: number of Hydrogens won't decrease by 1
                    esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + multiple);
                } else {
                    // in case of OH: number of Hydrogens won't decrease by 1
                    result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + multiple);
                }
                break;
            case "mercapto":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    // in case of SH: number of Hydrogens won't decrease by 1
                    esterSuffixResult.put("S", (esterSuffixResult.get("S") != null ? esterSuffixResult.get("S") : 0) + multiple);
                } else {
                    // in case of SH: number of Hydrogens won't decrease by 1
                    result.put("S", (result.get("S") != null ? result.get("S") : 0) + multiple);
                }
                break;
            case "imino":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put(NITROGEN, (esterSuffixResult.get(NITROGEN) != null ? esterSuffixResult.get(NITROGEN) : 0) + multiple);
                    // in case of NH: a Carbon loses 2 Hydrogens and gains 1, so it'll lose 1 Hydrogen in the end
                    esterSuffixNumberOfFunctionBounds += multiple;
                } else {
                    result.put(NITROGEN, (result.get(NITROGEN) != null ? result.get(NITROGEN) : 0) + multiple);
                    // in case of NH: a Carbon loses 2 Hydrogens and gains 1, so it'll lose 1 Hydrogen in the end
                    numberOfFunctionBounds += multiple;
                }
                break;
            case "amino":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put(NITROGEN, (esterSuffixResult.get(NITROGEN) != null ? esterSuffixResult.get(NITROGEN) : 0) + multiple);
                    // in case of NH2: a Carbon loses 1 Hydrogen and gains 2, so it'll gain 1 Hydrogen in the end
                    esterSuffixNumberOfFunctionBounds -= multiple;
                } else {
                    result.put(NITROGEN, (result.get(NITROGEN) != null ? result.get(NITROGEN) : 0) + multiple);
                    // in case of NH2: a Carbon loses 1 Hydrogen and gains 2, so it'll gain 1 Hydrogen in the end
                    numberOfFunctionBounds -= multiple;
                }
                break;
            case "arsino":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    // the same as in NH2, it'll gain 1 Hydrogen in the end
                    esterSuffixResult.put("As", (esterSuffixResult.get("As") != null ? esterSuffixResult.get("As") : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds -= multiple;
                } else {
                    // the same as in NH2, it'll gain 1 Hydrogen in the end
                    result.put("As", (result.get("As") != null ? result.get("As") : 0) + multiple);
                    numberOfFunctionBounds -= multiple;
                }
                break;
            case "phosphino":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    // the same as in NH2, it'll gain 1 Hydrogen in the end
                    esterSuffixResult.put("P", (esterSuffixResult.get("P") != null ? esterSuffixResult.get("P") : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds -= multiple;
                } else {
                    // the same as in NH2, it'll gain 1 Hydrogen in the end
                    result.put("P", (result.get("P") != null ? result.get("P") : 0) + multiple);
                    numberOfFunctionBounds -= multiple;
                }
                break;
            case "oxo":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + multiple);
                    esterSuffixNumberOfFunctionBounds += 2 * multiple;
                } else {
                    result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + multiple);
                    numberOfFunctionBounds += 2 * multiple;
                }
                break;
            case "formyl":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + multiple);
                    esterSuffixNumberOfCarbons += multiple;
                    esterSuffixNumberOfFunctionBounds += 2 * multiple;
                } else {
                    result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + multiple);
                    numberOfCarbons += multiple;
                    numberOfFunctionBounds += 2 * multiple;
                }
                break;
            case "phenyl":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixNumberOfCarbons += 6 * multiple;
                    // Each C6Hx takes away the same number of Hydrogens as 3 AlkeneBounds and a Cyclo do altogether
                    esterSuffixNumberOfAlkeneBounds += 4 * multiple;
                } else {
                    numberOfCarbons += 6 * multiple;
                    // Each C6Hx takes away the same number of Hydrogens as 3 AlkeneBounds and a Cyclo do altogether
                    numberOfAlkeneBounds += 4 * multiple;
                }
                break;
            case "amido":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + multiple);
                    esterSuffixResult.put(NITROGEN, (esterSuffixResult.get(NITROGEN) != null ? esterSuffixResult.get(NITROGEN) : 0) + multiple);
                    // in case of =O-C-NH2 a Carbon loses 3 Hydrogens and gains 2, so it'll lose 1 Hydro in the end
                    esterSuffixNumberOfFunctionBounds += multiple;
                } else {
                    result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + multiple);
                    result.put(NITROGEN, (result.get(NITROGEN) != null ? result.get(NITROGEN) : 0) + multiple);
                    // in case of =O-C-NH2 a Carbon loses 3 Hydrogens and gains 2, so it'll lose 1 Hydro in the end
                    numberOfFunctionBounds += multiple;
                }
                break;
            case "carboxy":
                if (esterSuffixMultiplier > 1 && esterSuffix) {
                    esterSuffixNumberOfCarbons += multiple;
                    esterSuffixResult.put(OXYGEN, (esterSuffixResult.get(OXYGEN) != null ? esterSuffixResult.get(OXYGEN) : 0) + 2 * multiple);
                    // each OH-C=O takes away 2 Hydrogens
                    esterSuffixNumberOfFunctionBounds += 2 * multiple;
                } else {
                    numberOfCarbons += multiple;
                    result.put(OXYGEN, (result.get(OXYGEN) != null ? result.get(OXYGEN) : 0) + 2 * multiple);
                    // each OH-C=O takes away 2 Hydrogens
                    numberOfFunctionBounds += 2 * multiple;
                }
                break;
            default:
        }
    }

    private class RelativeIndexHolder {
        public RelativeIndexHolder(int relIdx) {
            this.relIdx = relIdx;
        }

        private int relIdx;

        public int getRelIdx() {
            return relIdx;
        }

        public void setRelIdx(int relIdx) {
            this.relIdx = relIdx;
        }

        public void updateInsideBrackets(int startIdx, int endIdx) {
            relIdx -= startIdx;
        }

        public void updateOnShortenedString(int startIdx, int endIdx) {
            //------|---[---------]---|----
            if (relIdx >= startIdx && relIdx <= endIdx) {
                relIdx = startIdx;
            } else if (relIdx > endIdx) {
                relIdx -= (endIdx - startIdx + 1);
            }            
        } 
    }

    public static void main(String[] args) {
        //new ParseHer("2-1,1-dibromo-4-chloropentanylhexane"), // C: 6+5 or C: 11, H: 24, Br: 2, Cl: 1 // TODO: recursive at ramification  
        ParseHer parseHer = new ParseHer("1-[hexyl]arsino-9,17,17,17-tetrachloro-10-decyl-4-fluoro-9-heptyldihydroxy-2,11,12,14,15-pentaimino-13-iodo-13-nonadecyl-4-octyl-3,3-dioctoxyheptadecan-5,6,7,8,16-pentaone");
        System.out.println(String.format("Expression of %s is:", parseHer.molec));
        parseHer.parse().entrySet().stream()
            .sorted((entry1, entry2) -> entry1.getKey().compareTo(entry2.getKey()))
            .forEach(entry -> System.out.print(String.format("%s: %s, ", entry.getKey(), entry.getValue())));
    }
}

