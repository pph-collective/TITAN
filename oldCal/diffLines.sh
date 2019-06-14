#!/bin/bash
$1 = old
$2 = new
diff --unchanged-line-format="" --new-line-format="%dn: %L" --old-line-format="%dn: %L" $old $new

