use nom::{
    branch::{alt, permutation},
    bytes::complete::take,
    multi::{length_data, length_value, many0},
    number::complete::{be_u16, be_u32},
    sequence::tuple,
    IResult,
};

/// TODO: Rewrite it to be a lot simpler
/// Maybe initially forgo the enums and structs and just parse the data with functions

// Each snapgene packet has a header (5 bytes, tag + length of packet) and then the data which can be of a few packet "types"

#[derive(Debug, PartialEq)]
struct SnapGene;

struct SnapgenePacket {
    packet_header: PacketHeader,
    packet_data: PacketTypes,
}

#[derive(Debug, PartialEq)]
enum PacketTypes {
    Cookie(CookiePacket),
    Sequence(SequencePacket),
    Notes(NotesPacket),
    // Feature(FeaturePacket),
    // Primer(PrimerPacket),
}

// A packet header is the first 5 bytes of a packet
#[derive(Debug, PartialEq)]
struct PacketHeader {
    packet_tag: PacketType,
    packet_length: u32,
}

// Cookie is the initial packet
#[derive(Debug, PartialEq)]
struct CookiePacket {
    packet_header: PacketHeader,
    magic_cookie: [u8; 8],
    seq_type: SeqType,
    export_version: u16,
    import_version: u16,
}

// Appears once
#[derive(Debug, PartialEq)]
struct SequencePacket {
    packet_header: PacketHeader,
    sequence_flag: SequenceTopology,
    seq: Vec<u8>,
}

#[derive(Debug, PartialEq)]
struct NotesPacket {
    packet_header: PacketHeader,
    data: Vec<u8>,
}

impl NotesPacket {
    fn parse_notes(input: &[u8]) -> IResult<&[u8], Vec<u8>> {
        // Parse the packet header
        let (input, packet_header) = PacketHeader::parse_header(input)?;
        let (input, data) = take(packet_header.packet_length)(input)?;

        let combined = [data].concat();

        // Set the values
        let snap_parser = NotesPacket {
            packet_header,
            data: data.to_vec(),
        };
        Ok((input, combined))
    }
}
#[derive(Debug, PartialEq)]
enum SequenceTopology {
    Linear = 0x00, // Absent flag = linear
    Circular,      // Otherwise circular
}

#[derive(Debug, PartialEq)]
enum SeqType {
    DNA = 1,
    RNA = 2,
    Protein = 3,
}

#[derive(Debug, PartialEq)]
enum PacketType {
    Cookie = 0x09,
    Sequence = 0x00,
    Notes = 0x06,
    Feature = 0x0A,
    Primer = 0x0B,
    Other,
}

impl PacketHeader {
    fn parse_header(input: &[u8]) -> IResult<&[u8], PacketHeader> {
        // Header is tag (1 byte) + length (4 bytes) in big endian (long)
        let (input, (packet_tag, packet_length)) = tuple((take(1usize), be_u32))(input)?;

        let packet_header = PacketHeader {
            packet_tag: match packet_tag[0] {
                0x00 => PacketType::Sequence,
                0x06 => PacketType::Notes,
                0x09 => PacketType::Cookie,
                0x0A => PacketType::Feature,
                0x0B => PacketType::Primer,
                _ => PacketType::Other,
                // _ => panic!("Unknown packet type: byte = {:?}", packet_tag[0]),
            },
            packet_length,
        };
        // If packet type is other, remove the number of bytes from input
        if let PacketType::Other = packet_header.packet_tag {
            let (input, _) = take(packet_header.packet_length)(input)?;
            return Ok((input, packet_header));
        }

        Ok((input, packet_header))
    }
}

// This is a binary parser
impl CookiePacket {
    fn parse_cookie(input: &[u8]) -> IResult<&[u8], CookiePacket> {
        // Cookie always starts with 0x09, return otherwise
        if input[0] != 0x09 {
            panic!("Incorrect file or corrupted snapgene file");
        }

        // Parse the packet header (first 5 bytes), and then the rest of the cookie
        let (input, packet_header) = PacketHeader::parse_header(input)?;
        let (input, magic_cookie) = take(8usize)(input)?;
        let (input, seq_type) = be_u16(input)?; // (or take(2usize) and then convert to u16 manually)
        let (input, export_version) = be_u16(input)?;
        let (input, import_version) = be_u16(input)?;

        // Set the values
        let cookie_packet = CookiePacket {
            packet_header,
            magic_cookie: magic_cookie.try_into().unwrap(),
            seq_type: match seq_type {
                1 => SeqType::DNA,
                2 => SeqType::RNA,
                3 => SeqType::Protein,
                _ => panic!("Unknown sequence type"),
            },
            export_version,
            import_version,
        };

        Ok((input, cookie_packet))
    }
}

impl SequencePacket {
    fn parse_dnapacket(input: &[u8]) -> IResult<&[u8], Vec<u8>> {
        // Parse the packet header
        let (input, packet_header) = PacketHeader::parse_header(input)?;
        let (input, sequence_flag) = take(1usize)(input)?;

        // Sequence length is packet length - 1 (due to sequence flag)
        let length = packet_header.packet_length - 1;
        let (input, seq) = take(length as usize)(input)?;
        let combined = [sequence_flag, seq].concat();

        // Set the values
        let snap_parser = SequencePacket {
            packet_header,
            sequence_flag: match sequence_flag[0] {
                0x00 => SequenceTopology::Linear,
                _ => SequenceTopology::Circular,
            },
            seq: String::from_utf8(seq.to_vec()).unwrap().into_bytes(),
        };
        Ok((input, combined))
    }
}

impl SnapGene {
    fn parse(input: &[u8]) -> IResult<&[u8], SnapGene> {
        let (input, cookie) = CookiePacket::parse_cookie(input).unwrap();

        // The rest of the inputs can be the any of the PacketTypes in any order, so we need to parse them all
        let (input, packets) = many0(alt((
            SequencePacket::parse_dnapacket,
            NotesPacket::parse_notes,
        )))(input)
        .unwrap();

        println!("{:?}", packets);

        Ok((input, SnapGene))
    }
}

fn main() {
    let f = std::fs::read("./pcDNA4 TO.dna").unwrap();

    // Cookie packet is always first
    // Parse the rest of the packets
    let (input, snapgene) = SnapGene::parse(&f).unwrap();
}
